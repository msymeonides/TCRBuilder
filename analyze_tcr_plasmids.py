from __future__ import annotations
import argparse
import csv
import re
import statistics
import sys
import zipfile
from collections import defaultdict, Counter
from pathlib import Path
from xml.etree import ElementTree as ET

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

# ====================== USER CONFIGURATION =============================
PURE_MEDIAN_THRESHOLD = 98.0     # strictly greater than 98% is 'pure'
BACKBONE_WINDOW_BP = 500
BACKBONE_MIN_MEDIAN_CONSENSUS = 99.0  # flag windows STRICTLY below 99%
USABLE_MEDIAN_THRESHOLD = 90.0   # strictly greater than 90% can be purified
REQUIRED_IDENTITY_PERCENT = 100.0 # strict nucleotide identity
MIN_EXPECTED_COVERAGE = 100.0   # require full-length intended eBlock
# The ONLY exception: TCRb eBlock junction, 5th base from 3' end, G -> A.
ALLOW_TCRB_JUNCTION_G_TO_A = True
ALLOW_SYNONYMOUS_SUBSTITUTIONS = True
TRANSLATION_START_MOTIF = "GGGGATCGCTCGAGACCATG"  # ATG is the last 3 bases
# Diagnostics, NOT permission to classify an incorrect sequence as usable.
DIAGNOSTIC_MIN_IDENTITY = 95.0
DIAGNOSTIC_KMER_SIZE = 15
DIAGNOSTIC_TOP_CANDIDATES = 8
# Explicitly failed sequencing is treated as no plasmid, per user assumption.
FAIL_MEANS_NO_PLASMID = True
# Each sequencing submission/plate contains ONE chain. Add new IDs here.
# IDs are extracted from the DATE_PLATEID dataset folder name.
PLATE_CHAIN_MAP = {
    '1198368': 'TCRa',
    '1198373': 'TCRb',
}
# Local alignment penalties: prefer a nearly complete matching eBlock.
ALIGN_MATCH, ALIGN_MISMATCH, ALIGN_GAP_OPEN, ALIGN_GAP_EXTEND = 2, -3, -6, -1
# =========================================================================

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
      'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
RELOCATIONS = {
    'D1':'A9','D2':'A10','D3':'B9','D4':'C9','D5':'C10','D6':'C11',
    'H1':'E9','H2':'E10','H3':'F9','H4':'G9','H5':'G10','H6':'G11'
}
COLONY_FIELDS = ['data_folder','plate','sample_number','sample_name','colony','assembly_well',
 'idt_well','insert_name','sequencing_status','sequencing_size','fasta_file','per_base_file','status','reason',
 'median_consensus_pct','mean_consensus_pct','minimum_consensus_pct',
 'bases_below_98pct','bases_below_90pct','region_start_1based','region_end_1based',
 'expected_length_bp','aligned_expected_bp','coverage_pct','identity_pct',
 'matches','mismatches','insertions_bp','deletions_bp','strand',
 'best_other_insert','best_other_identity_pct','best_other_coverage_pct','allowed_junction_variant','synonymous_mismatches','synonymous_changes']
INSERT_FIELDS = ['data_folder','chain','idt_well','assembly_well','insert_name',
 'status','preferred_sample','preferred_colony','preferred_median_consensus_pct',
 'preferred_identity_pct','colonies_on_layout','colonies_with_files',
 'verified_colonies','purify_colonies','other_or_failed_colonies','sequencing_failed_colonies','notes']


def xlsx_rows(path):
    """Read simple workbook sheets without Excel-specific third-party dependencies."""
    with zipfile.ZipFile(path) as z:
        strings = []
        if 'xl/sharedStrings.xml' in z.namelist():
            root = ET.fromstring(z.read('xl/sharedStrings.xml'))
            strings = [''.join(t.text or '' for t in si.iter('{%s}t' % NS['m']))
                       for si in root]
        sheets = sorted((p for p in z.namelist() if re.fullmatch(r'xl/worksheets/sheet\d+\.xml', p)),
                        key=lambda s: int(re.search(r'sheet(\d+)',s).group(1)))
        for path_in_zip in sheets:
            root = ET.fromstring(z.read(path_in_zip))
            rows = []
            for row in root.findall('.//m:sheetData/m:row', NS):
                cells = {}
                for c in row.findall('m:c', NS):
                    column = re.match(r'[A-Z]+',c.attrib['r']).group()
                    idx = 0
                    for letter in column:idx = idx * 26 + ord(letter) - 64
                    v = c.find('m:v',NS)
                    if c.attrib.get('t') == 'inlineStr':
                        is_el = c.find('m:is',NS)
                        value = ''.join(t.text or '' for t in is_el.iter('{%s}t' % NS['m'])) if is_el is not None else ''
                    elif v is None: value = ''
                    elif c.attrib.get('t') == 's':value = strings[int(v.text)]
                    else:value = v.text or ''
                    cells[idx-1] = value
                if cells:
                    output = [''] * (max(cells)+1)
                    for idx,value in cells.items():output[idx] = value
                    rows.append(output)
            yield rows


def build_layouts(layout_dir):
    idt_file = layout_dir/'IDT-eBlocks_v4.xlsx'
    plate_file = layout_dir/'Miniprep sequencing.xlsx'
    if not idt_file.exists() or not plate_file.exists():
        raise FileNotFoundError('Place both IDT-eBlocks_v4.xlsx and Miniprep sequencing.xlsx in layouts/')
    eblocks = {}
    for sheet in xlsx_rows(idt_file):
        for r in sheet[1:]:
            if len(r)>=3 and re.fullmatch(r'[A-H](?:[1-9]|1[0-2])',r[0].strip().upper()) and r[2].strip():
                well, name, seq = r[:3]
                seq = re.sub(r'\s+','',seq).upper()
                if not re.fullmatch('[ACGT]+',seq): raise ValueError('Invalid eBlock DNA: '+name)
                eblocks[well.strip().upper()] = (name.strip(),seq)
    layouts = []
    # Only the first two worksheets are sequencing sample lists.
    # Worksheets 3 and 4 describe physical miniprep well positions and
    # must NOT be interpreted as additional sequencing plates.
    for sheet_idx, sheet in enumerate(xlsx_rows(plate_file)):
        if sheet_idx >= 2:
            break
        plate = 'TCRa' if sheet_idx == 0 else 'TCRb'
        for r in sheet[1:]:
            if len(r)<2:continue
            number = str(r[0]).strip(); name=str(r[1]).strip().upper()
            m = re.fullmatch(r'([A-H](?:1[0-2]|[1-9]))C(\d+)',name)
            if not m:continue
            well,colony=m.group(1),m.group(2)
            original = RELOCATIONS.get(well,well)
            if original not in eblocks:raise ValueError(f'{plate} {name}: unmapped original {original}')
            insert_name,insert_seq = eblocks[original]
            chain = 'TCRa' if original[0] in 'ABCD' else 'TCRb'
            if chain != plate:raise ValueError(f'Unexpected chain for {plate} {name}')
            layouts.append(dict(plate=plate,sample_number=int(float(number)),sample_name=name,
                colony=int(colony),assembly_well=well,idt_well=original,
                insert_name=insert_name,expected_seq=insert_seq))
    if len(eblocks)!=60:print(f'WARNING: expected 60 IDT sequences, found {len(eblocks)}',file=sys.stderr)
    return eblocks, layouts


def get_seq(path):
    lines=path.read_text(encoding='utf-8-sig').splitlines()
    seq=''.join(line.strip() for line in lines if line.strip() and not line.lstrip().startswith('>')).upper()
    seq=re.sub(r'\s+','',seq)
    if not seq or not re.fullmatch('[ACGTN]+',seq):raise ValueError(f'Unrecognized FASTA: {path}')
    return seq


def get_bases(path):
    positions=[]; seq=[]; percentages=[]
    with path.open(newline='',encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            try:
                pos=int(row['Position']);base=row['Consensus Base'].strip().upper()
                percent=float(row['Consensus (%)'])
            except (ValueError,KeyError,TypeError):continue
            if base not in 'ACGTN' or len(base)!=1:continue
            positions.append(pos);seq.append(base);percentages.append(percent)
    if not seq:raise ValueError('No usable base calls in '+str(path))
    return ''.join(seq),positions,percentages


def make_aligner():
    aligner=PairwiseAligner(mode='local')
    aligner.match_score=ALIGN_MATCH
    aligner.mismatch_score=ALIGN_MISMATCH
    aligner.open_gap_score=ALIGN_GAP_OPEN
    aligner.extend_gap_score=ALIGN_GAP_EXTEND
    return aligner


def evaluate(aligner, reference, plasmid, include_query_indices=False):
    """Align intended eBlock (reference) to entire consensus circular plasmid.
    Duplicate the first ~reference length bases to allow origin-spanning alignments.
    Return reference coverage, mismatches, indels, and subject indexes.
    """
    n=len(plasmid)
    subject=plasmid+plasmid[:min(n,len(reference)+30)]
    options=[]
    for strand,query in [('+',reference),('-',str(Seq(reference).reverse_complement()))]:
        al=aligner.align(subject,query)[0]
        pairs=al.indices
        matches=mismatches=insertions=deletions=0
        aligned_ref=set();mapped=[]
        subject_positions=[]
        for s,q in zip(pairs[0],pairs[1]):
            if q>=0:aligned_ref.add(int(q))
            if s>=0 and q>=0:
                subject_positions.append(int(s)%n)
                mapped.append(int(s)%n)
                if subject[int(s)]==query[int(q)]:matches+=1
                else:mismatches+=1
            elif q>=0:deletions+=1
            elif s>=0:insertions+=1
        coverage=100*len(aligned_ref)/len(reference)
        identity=100*matches/max(1,matches+mismatches+insertions+deletions)
        options.append(dict(strand=strand,coverage_pct=coverage,identity_pct=identity,
            matches=matches,mismatches=mismatches,insertions_bp=insertions,deletions_bp=deletions,
            aligned_expected_bp=len(aligned_ref),query_indices=sorted(set(mapped)),
            region_start_1based=min(subject_positions)+1 if subject_positions else None,
            region_end_1based=max(subject_positions)+1 if subject_positions else None))
    return max(options,key=lambda d:(d['coverage_pct']>=MIN_EXPECTED_COVERAGE,d['identity_pct'],d['coverage_pct'],d['matches']))


def accepted_reference_sequences(well, reference):
    """Permit only the specific TCRb G->A substitution at index len(ref)-5."""
    variants = [reference]
    if ALLOW_TCRB_JUNCTION_G_TO_A and well[0] in 'EFGH':
        idx = len(reference) - 5
        if idx < 0 or reference[idx] != 'G':
            raise ValueError(f'TCRb junction exception cannot apply to {well}')
        variants.append(reference[:idx] + 'A' + reference[idx + 1:])
    return variants


def translation_start(reference):
    """Find the unique motif; last ATG starts the ORF (0-based index)."""
    positions=[m.start() for m in re.finditer(TRANSLATION_START_MOTIF,reference)]
    if len(positions)!=1:
        raise ValueError(f'Expected one ATG anchor motif, found {len(positions)}')
    start=positions[0]+len(TRANSLATION_START_MOTIF)-3
    if reference[start:start+3]!='ATG':raise ValueError('Translation motif ATG mismatch')
    return start


def classify_reference_match(well,reference,plasmid,aligner=None):
    """Accept full-length ungapped matches with only synonymous coding substitutions
    and/or the specifically permitted TCRb G->A homology-arm difference.
    Return actual observed reference-oriented insert sequence and auditable edits.
    """
    if aligner is None:aligner=make_aligner()
    n=len(plasmid)
    if not n:return None
    subject=plasmid+plasmid[:min(n,len(reference)+30)]
    starts=translation_start(reference)
    best=None
    for strand,query in (('+',reference),('-',str(Seq(reference).reverse_complement()))):
        alignment=aligner.align(subject,query)[0]
        indices=alignment.indices
        # Require ALL reference bases in one contiguous ungapped alignment.
        pairs=[(int(u),int(v)) for u,v in zip(indices[0],indices[1])]
        if len(pairs)!=len(reference) or any(u<0 or v<0 for u,v in pairs):continue
        if sorted(v for _,v in pairs)!=list(range(len(reference))):continue
        if any(pairs[i+1][0]!=pairs[i][0]+1 or pairs[i+1][1]!=pairs[i][1]+1 for i in range(len(pairs)-1)):continue
        observed=''.join(subject[u] for u,v in pairs)
        if strand=='-':observed=str(Seq(observed).reverse_complement())
        if any(base not in 'ACGT' for base in observed):continue
        differences=[i for i,(a,b) in enumerate(zip(reference,observed)) if a!=b]
        if not differences:
            return dict(sequence=observed,allowed_junction_variant=False,
                        synonymous_mismatches=0,synonymous_changes='',strand=strand)
        junction= len(reference)-5
        allowed_junction=[];coding_differences=[]
        for i in differences:
            if (ALLOW_TCRB_JUNCTION_G_TO_A and well[0] in 'EFGH'
                and i==junction and reference[i]=='G' and observed[i]=='A'):
                allowed_junction.append(i)
            else:
                coding_differences.append(i)
        if coding_differences and not ALLOW_SYNONYMOUS_SUBSTITUTIONS:continue
        if any(i<starts or i+ (3-((i-starts)%3))>len(reference) for i in coding_differences):continue
        synonymous=[];valid=True
        for codon_start in sorted({starts+3*((i-starts)//3) for i in coding_differences}):
            ref_codon=reference[codon_start:codon_start+3]
            obs_codon=observed[codon_start:codon_start+3]
            if len(ref_codon)!=3 or len(obs_codon)!=3 or str(Seq(ref_codon).translate()) != str(Seq(obs_codon).translate()):
                valid=False;break
        if not valid:continue
        for i in coding_differences:
            codon_start=starts+3*((i-starts)//3)
            synonymous.append(f'{reference[i]}{i+1}{observed[i]} '
                              f'(codon {((i-starts)//3)+1}: '
                              f'{reference[codon_start:codon_start+3]}>'
                              f'{observed[codon_start:codon_start+3]})')
        candidate=dict(sequence=observed,allowed_junction_variant=bool(allowed_junction),
                       synonymous_mismatches=len(coding_differences),
                       synonymous_changes='; '.join(synonymous),strand=strand)
        if best is None or len(coding_differences)<best['synonymous_mismatches']:
            best=candidate
    return best


def matching_reference(well,reference,plasmid,aligner=None):
    """Return observed insert if it is an approved translation-equivalent variant."""
    result=classify_reference_match(well,reference,plasmid,aligner)
    return result['sequence'] if result else None


def locate_folder(data_dir, phrase):
    exact=data_dir/phrase
    if exact.is_dir():return exact
    found=[p for p in data_dir.rglob('*') if p.is_dir() and phrase.lower() in p.name.lower()]
    if len(found)!=1:raise FileNotFoundError(f'Expected one {phrase} folder beneath {data_dir}; found {len(found)}')
    return found[0]


def find_sample_file(folder,sample_number,sample_name,extension):
    files=list(folder.rglob('*'+extension))
    # Sequencing provider often prefixes zero-padded sample number, then identifier.
    starts=[p for p in files if re.match(rf'^0*{sample_number}(?!\d)',p.stem,re.I)]
    named=[p for p in starts if sample_name.lower() in p.stem.lower()]
    # Never assign another colony's file just because its plate-local number matches.
    # Number-only filenames are acceptable only if truly unlabelled.
    other_named=[p for p in starts if re.search(r'[A-H](?:1[0-2]|[1-9])C\d+',p.stem,re.I)]
    candidates=named or ([] if other_named else starts)
    if len(candidates)==1:return candidates[0]
    if len(candidates)>1:raise ValueError(f'Ambiguous {extension} for {sample_number} {sample_name}: {candidates}')
    return None


def output_csv(path,fields,records):
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fields,extrasaction='ignore');w.writeheader();w.writerows(records)


def plate_id_from_folder(data_dir):
    match = re.fullmatch(r'\d{8}_([A-Za-z0-9_-]+)', data_dir.name)
    if not match:
        raise ValueError(f'Data directory must follow DATE_PLATEID format: {data_dir.name}')
    return match.group(1)


def select_chain(data_dir, chain_map):
    plate_id = plate_id_from_folder(data_dir)
    if plate_id not in chain_map:
        raise ValueError(
            f'Unknown sequencing plate ID {plate_id!r} ({data_dir.name}). '
            'Add it to PLATE_CHAIN_MAP at the top of the script or pass '
            f'--plate-chain {plate_id}=TCRa (or =TCRb). '
            'No chain was inferred because that could produce incorrect classifications.'
        )
    chain = chain_map[plate_id]
    if chain not in ('TCRa', 'TCRb'):
        raise ValueError(f'Invalid chain {chain!r} for plate ID {plate_id}')
    return plate_id, chain


def load_sequencing_summary(data_dir, plate_id):
    """Load the provider's per-plate sequencing success/failure manifest."""
    path = data_dir / f'{plate_id}_summary.csv'
    if not path.is_file():
        print(f'  WARNING: {path.name} absent; sequencing failures cannot be distinguished from missing data', file=sys.stderr)
        return None
    result = {}
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        required = {'Sample #', 'Sample Name', 'Sample Status'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'{path}: required columns missing: {sorted(required-set(reader.fieldnames or []))}')
        for record in reader:
            num = str(record['Sample #']).strip()
            name = str(record['Sample Name']).strip().upper()
            if not num.isdigit() or not name:
                raise ValueError(f'{path}: invalid sample identifier: {num!r} {name!r}')
            key = int(num)
            if key in result:
                raise ValueError(f'{path}: duplicate sample number {key}')
            result[key] = dict(name=name, status=record['Sample Status'].strip(), size=(record.get('Sample Size') or '').strip())
    return result


def process_folder(data_dir,eblocks,layouts,chain_map):
    plate_id, chain = select_chain(data_dir, chain_map)
    layouts = [sample for sample in layouts if sample['plate'] == chain]
    eblocks = {well: info for well, info in eblocks.items()
               if ('TCRa' if well[0] in 'ABCD' else 'TCRb') == chain}
    print(f'  Plate {plate_id} -> {chain}: {len(layouts)} colony samples, '
          f'{len(eblocks)} expected inserts', flush=True)
    per_folder=locate_folder(data_dir,'PER_BASE_DATA')
    fasta_folder=locate_folder(data_dir,'FASTA_FILES')
    output_dir=data_dir/'results';output_dir.mkdir(exist_ok=True)
    sequencing_summary=load_sequencing_summary(data_dir,plate_id)
    aligner=make_aligner()
    colony_results=[]
    for ix,sample in enumerate(layouts,1):
        row={k:v for k,v in sample.items() if k!='expected_seq'}
        row['data_folder']=data_dir.name
        manifest_entry = sequencing_summary.get(sample['sample_number']) if sequencing_summary is not None else None
        if sequencing_summary is not None:
            if manifest_entry is None:
                row['sequencing_status']='not in summary'
                row['status']='missing sequencing summary entry'
                row['reason']='Sample is in layout but absent from provider summary; review mapping/submission'
                colony_results.append(row)
                continue
            if manifest_entry['name'] != sample['sample_name']:
                row['sequencing_status']=manifest_entry['status']
                row['sequencing_size']=manifest_entry['size']
                row['status']='sequencing summary mismatch'
                row['reason']=f"Summary sample name {manifest_entry['name']} differs from layout {sample['sample_name']}"
                colony_results.append(row)
                continue
            row['sequencing_status']=manifest_entry['status']
            row['sequencing_size']=manifest_entry['size']
            if manifest_entry['status'].casefold() == 'fail':
                row['status']='sequencing failed'
                row['reason']='Provider reported Fail; resequence this colony if retained'
                colony_results.append(row)
                continue
            if manifest_entry['status'].casefold() != 'success':
                row['status']='unknown sequencing status'
                row['reason']=f"Unrecognized provider status: {manifest_entry['status']}"
                colony_results.append(row)
                continue
        else:
            row['sequencing_status']='summary unavailable'
        fasta=find_sample_file(fasta_folder,sample['sample_number'],sample['sample_name'],'.fasta')
        if fasta is None:fasta=find_sample_file(fasta_folder,sample['sample_number'],sample['sample_name'],'.fa')
        basefile=find_sample_file(per_folder,sample['sample_number'],sample['sample_name'],'.csv')
        row['fasta_file']=str(fasta.relative_to(data_dir)) if fasta else ''
        row['per_base_file']=str(basefile.relative_to(data_dir)) if basefile else ''
        row['expected_length_bp']=len(sample['expected_seq'])
        if not fasta or not basefile:
            row['status']='missing data';row['reason']='Missing '+('FASTA ' if not fasta else '')+('per-base CSV' if not basefile else '')
            colony_results.append(row);continue
        try:
            fasta_seq=get_seq(fasta)
            consensus,positions,percentages=get_bases(basefile)
            # Locate purity interval using per-base consensus data, never by assuming plasmid offset.
            approved=classify_reference_match(sample['idt_well'],sample['expected_seq'],fasta_seq,aligner)
            accepted=approved['sequence'] if approved else None
            comparison_reference=accepted or sample['expected_seq']
            row['allowed_junction_variant']=bool(approved and approved['allowed_junction_variant'])
            row['synonymous_mismatches']=approved['synonymous_mismatches'] if approved else 0
            row['synonymous_changes']=approved['synonymous_changes'] if approved else ''
            base_alignment=evaluate(aligner,comparison_reference,consensus)
            fasta_alignment=evaluate(aligner,comparison_reference,fasta_seq)
            indices=base_alignment['query_indices']
            chosen=[percentages[j] for j in indices if j<len(percentages)]
            if not chosen:raise ValueError('Could not locate candidate insert in per-base data')
            row.update({k:v for k,v in fasta_alignment.items() if k not in ('query_indices','region_start_1based','region_end_1based')})
            row.update(median_consensus_pct=round(statistics.median(chosen),4),
                mean_consensus_pct=round(statistics.mean(chosen),4),
                minimum_consensus_pct=min(chosen),bases_below_98pct=sum(x<98 for x in chosen),
                bases_below_90pct=sum(x<90 for x in chosen),
                region_start_1based=min(positions[j] for j in indices),
                region_end_1based=max(positions[j] for j in indices))
            complete=fasta_alignment['coverage_pct'] >= MIN_EXPECTED_COVERAGE - 1e-7
            matched=bool(approved) and complete and fasta_alignment['identity_pct'] >= REQUIRED_IDENTITY_PERCENT - 1e-7
            # Purity threshold cannot be applied to a short/truncated alignment.
            located=base_alignment['coverage_pct'] >= MIN_EXPECTED_COVERAGE - 1e-7
            pure=row['median_consensus_pct']>PURE_MEDIAN_THRESHOLD
            usable=row['median_consensus_pct']>USABLE_MEDIAN_THRESHOLD
            if matched and located and pure:
                row['status']='verified good and pure';row['reason']='Full-length accepted insert; synonymous substitutions: '+str(row['synonymous_changes'])
            elif matched and located and usable:
                row['status']='good but need to purify';row['reason']='Full-length accepted insert; synonymous substitutions: '+str(row['synonymous_changes'])
            else:
                row['status']='inadequate'
                row['reason']='; '.join(x for x in [
                    'Incomplete FASTA alignment' if not complete else '',
                    'FASTA insert does not exactly match expected' if not matched else '',
                    'Incomplete per-base insert alignment' if not located else '',
                    'Median consensus at/below usable cutoff' if not usable else ''] if x)
            if not matched:
                # Report definite alternative insert matches without thousands of costly alignments.
                # A 100% full-length exact match to another eBlock strongly suggests a mix-up.
                other_matches=[]
                for well,(name,seq) in eblocks.items():
                    if well==sample['idt_well']:continue
                    if matching_reference(well,seq,fasta_seq):
                        other_matches.append(name)
                if other_matches:
                    row['best_other_insert']='; '.join(other_matches)
                    row['best_other_identity_pct']=100.0
                    row['best_other_coverage_pct']=100.0

        except Exception as exc:
            row['status']='analysis error';row['reason']=f'{type(exc).__name__}: {exc}'
        colony_results.append(row)
        if ix%25==0:print(f'  analyzed {ix}/{len(layouts)} colonies',flush=True)
    by_well=defaultdict(list)
    for r in colony_results:by_well[r['idt_well']].append(r)
    insert_results=[]
    priority={'verified good and pure':2,'good but need to purify':1}
    for well,(name,seq) in eblocks.items():
        samples=by_well.get(well,[])
        good=[r for r in samples if r['status'] in priority]
        good.sort(key=lambda r:(priority[r['status']],r.get('median_consensus_pct') or 0,r.get('sample_number') or 0),reverse=True)
        preferred=good[0] if good else None
        if preferred:status=preferred['status']
        elif not samples:status='need to retransform - no sequences'
        elif any(r['status'] in ('missing data','analysis error','missing sequencing summary entry',
                                 'sequencing summary mismatch','unknown sequencing status') for r in samples):
            status='incomplete sequencing data - review manually'
        elif any(r['status']=='sequencing failed' for r in samples):
            status='need to resequence - failed colonies'
        elif all(r['status']=='missing data' for r in samples):
            status='incomplete sequencing data - review manually'
        else:status='need to retransform - sequences all inadequate'
        row=dict(data_folder=data_dir.name,chain='TCRa' if well[0] in 'ABCD' else 'TCRb',
            idt_well=well,assembly_well=next((k for k,v in RELOCATIONS.items() if v==well),well),
            insert_name=name,status=status,preferred_sample=preferred['sample_name'] if preferred else '',
            preferred_colony=preferred['colony'] if preferred else '',
            preferred_median_consensus_pct=preferred.get('median_consensus_pct','') if preferred else '',
            preferred_identity_pct=preferred.get('identity_pct','') if preferred else '',
            colonies_on_layout=len(samples),
            colonies_with_files=sum(bool(r.get('fasta_file')) and bool(r.get('per_base_file')) for r in samples),
            verified_colonies='; '.join(r['sample_name'] for r in samples if r['status']=='verified good and pure'),
            purify_colonies='; '.join(r['sample_name'] for r in samples if r['status']=='good but need to purify'),
            other_or_failed_colonies='; '.join(r['sample_name'] for r in samples if r['status'] not in priority),
            sequencing_failed_colonies='; '.join(r['sample_name'] for r in samples if r['status']=='sequencing failed'),
            notes='; '.join(r['sample_name']+': '+r.get('reason','') for r in samples if r['status'] in ('missing data','analysis error','sequencing failed','missing sequencing summary entry','sequencing summary mismatch','unknown sequencing status')))
        insert_results.append(row)
    output_csv(output_dir/'colony_results.csv',COLONY_FIELDS,colony_results)
    output_csv(output_dir/'insert_summary.csv',INSERT_FIELDS,insert_results)
    output_csv(output_dir/'verified_good_and_pure.csv',INSERT_FIELDS,
               [r for r in insert_results if r['status']=='verified good and pure'])
    output_csv(output_dir/'action_required.csv',INSERT_FIELDS,
               [r for r in insert_results if r['status']!='verified good and pure'])
    counts=defaultdict(int)
    for row in insert_results:counts[row['status']]+=1
    with (output_dir/'summary.txt').open('w',encoding='utf-8') as f:
        f.write(f'Dataset: {data_dir.name}\nPlate ID: {plate_id}\nChain: {chain}\nColonies on layout: {len(layouts)}\n'
                f'Intended inserts: {len(eblocks)}\n'
                f'Sequencing summary: {plate_id}_summary.csv ({"loaded" if sequencing_summary is not None else "not found"})\n'
                f'Provider sequencing failures: {sum(r["status"] == "sequencing failed" for r in colony_results)}\n'
                f'Pure > {PURE_MEDIAN_THRESHOLD}%, usable > {USABLE_MEDIAN_THRESHOLD}%, '
                f'identity >= {REQUIRED_IDENTITY_PERCENT}%, coverage >= {MIN_EXPECTED_COVERAGE}%\n\n')
        for status,n in sorted(counts.items()):f.write(f'{status}: {n}\n')
    print(f'  Finished {data_dir.name}: '+', '.join(f'{k}: {n}' for k,n in counts.items()))
    return output_dir


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--diagnostics',action='store_true',help='Also write one optional detailed colony CSV')
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent,
                        help='Working directory containing layouts/ and DATE_PLATEID/ folders')
    parser.add_argument('--data-folder',action='append',help='Optional folder name, can repeat')
    parser.add_argument('--plate-chain', action='append', default=[], metavar='PLATEID=CHAIN',
                        help='Map an additional/overridden sequencing plate ID to TCRa or TCRb; repeatable')
    a=parser.parse_args();root=a.root.resolve()
    chain_map = dict(PLATE_CHAIN_MAP)
    for item in a.plate_chain:
        if '=' not in item:
            parser.error('--plate-chain must be PLATEID=TCRa or PLATEID=TCRb')
        plate_id, chain = (part.strip() for part in item.split('=', 1))
        if not plate_id or chain not in ('TCRa', 'TCRb'):
            parser.error(f'Invalid --plate-chain {item!r}; expected PLATEID=TCRa or PLATEID=TCRb')
        chain_map[plate_id] = chain
    eblocks,layouts=build_layouts(root/'layouts')
    folders=[root/f for f in a.data_folder] if a.data_folder else sorted(
        p for p in root.iterdir() if p.is_dir() and re.fullmatch(r'\d{8}_[A-Za-z0-9_-]+',p.name))
    if not folders:parser.error('No DATE_PLATEID data folders found under '+str(root))
    for folder in folders:
        print('Processing '+str(folder),flush=True)
        write_actionable_reports(folder, eblocks, layouts, chain_map, diagnostics=a.diagnostics)



# =================== INTEGRATED DIAGNOSTIC AND ACTION REPORTS ==================
# This stage leaves the original files intact and writes the four requested action
# lists. It compares exact sequences against ALL 60 inserts; nearest-neighbor
# alignments are only diagnostic and never qualify a colony as ready.
ACTION_FIELDS = ['chain','idt_well','assembly_well','insert_name','action',
                 'selected_sample','sequencing_well','origin_assembly_well',
                 'median_consensus_pct','identity_pct','coverage_pct',
                 'other_candidate_samples','synonymous_mismatches','synonymous_changes','reason']
DETAIL_FIELDS = ['sample_number','sample_name','original_assembly_well','original_idt_well',
                 'sequencing_status','fasta_file','per_base_file','exact_reference_wells',
                 'exact_other_chain_wells','best_candidate_well','best_candidate_name',
                 'best_candidate_identity_pct','best_candidate_coverage_pct',
                 'best_candidate_mismatches','best_candidate_insertions_bp',
                 'best_candidate_deletions_bp','median_consensus_pct',
                 'per_base_coverage_pct','plasmid_integrity','integrity_flagged_windows','backbone_windows_evaluated','allowed_junction_variant','synonymous_mismatches','synonymous_changes','decision','notes']

def _exact_refs(seq, eblocks, aligner):
    # K-mer shortlist avoids an expensive all-versus-all whole-plasmid alignment.
    candidates=_rank_kmers(seq,eblocks)
    return [well for well in candidates if matching_reference(well,eblocks[well][1],seq,aligner) is not None]

def _rank_kmers(subject, eblocks):
    k=DIAGNOSTIC_KMER_SIZE
    subject=subject+subject[:600]
    kmers=set(subject[i:i+k] for i in range(max(0,len(subject)-k+1)))
    scored=[]
    for well,(_,ref) in eblocks.items():
        rev=str(Seq(ref).reverse_complement())
        for q in (ref,rev):
            score=sum(q[i:i+k] in kmers for i in range(max(0,len(q)-k+1)))
            scored.append((score,well))
    scores={}
    for val,well in scored:scores[well]=max(scores.get(well,0),val)
    return sorted(scores, key=lambda well:(scores[well],well),reverse=True)[:DIAGNOSTIC_TOP_CANDIDATES]

def _consensus_for_reference(aligner, ref, consensus, percentages):
    if consensus is None:return (None,None)
    a=evaluate(aligner,ref,consensus)
    if a['coverage_pct'] < MIN_EXPECTED_COVERAGE - 1e-7:return (None,a['coverage_pct'])
    vals=[percentages[i] for i in a['query_indices'] if i<len(percentages)]
    return (round(statistics.median(vals),4) if vals else None,a['coverage_pct'])

def assess_backbone_integrity(aligner, insert_reference, consensus, positions, percentages):
    """500-base windows entirely outside insert, in per-base positional order.

    Coordinate ranges are inclusive, 1-based provider Position values. Windows
    are created separately for each contiguous non-insert run. A short final
    window is still evaluated. The report is a consensus QC flag, NOT structural
    confirmation that the complete plasmid sequence is correct.
    """
    if not consensus or len(consensus) != len(percentages) or len(consensus) != len(positions):
        return 'not assessed', 'Incomplete per-base data', []
    mapping = evaluate(aligner, insert_reference, consensus)
    if (mapping['coverage_pct'] < MIN_EXPECTED_COVERAGE-1e-7
            or len(mapping['query_indices']) != len(insert_reference)):
        return 'not assessed', 'Insert boundary could not be located reliably', []
    masked = set(mapping['query_indices'])
    # If the insert spans the circular origin, masked positions naturally
    # occur on both ends. No excluded base is permitted into a window.
    regions=[];current=[]
    for idx in sorted(range(len(positions)), key=lambda j:positions[j]):
        if idx in masked:
            if current: regions.append(current);current=[]
            continue
        if current and positions[idx] != positions[current[-1]]+1:
            regions.append(current);current=[]
        current.append(idx)
    if current: regions.append(current)
    windows=[]
    for region in regions:
        for start in range(0,len(region),BACKBONE_WINDOW_BP):
            idxs=region[start:start+BACKBONE_WINDOW_BP]
            value=statistics.median(percentages[i] for i in idxs)
            windows.append((positions[idxs[0]], positions[idxs[-1]], round(value,3)))
    if not windows:return 'not assessed','No backbone bases outside insert',[]
    flagged=[w for w in windows if w[2]<BACKBONE_MIN_MEDIAN_CONSENSUS]
    details='; '.join(f'{a}-{b}: {v:.2f}%' for a,b,v in flagged)
    return ('follow up' if flagged else 'intact'),details,windows


def read_miniprep_mapping(layout_dir, chain, samples):
    """Physical colony -> miniprep plate coordinate; fail visibly on layout mismatch."""
    path=layout_dir/'Miniprep sequencing.xlsx'
    if not path.exists(): return {}, 'Miniprep layout workbook is missing.'
    sheets=list(xlsx_rows(path))
    sheet_num=2 if chain=='TCRa' else 3
    if len(sheets)<=sheet_num:return {}, f'{chain} miniprep worksheet is missing.'
    matrix=sheets[sheet_num]
    mapping={}
    for row_idx,row in enumerate(matrix[1:9]):
        for col_idx,raw in enumerate(row[1:13]):
            colony=str(raw).strip().upper()
            if not re.fullmatch(r'[A-H](?:1[0-2]|[1-9])C\d+',colony):continue
            if colony in mapping: return {}, f'Duplicate colony in {chain} miniprep sheet: {colony}.'
            mapping[colony]='ABCDEFGH'[row_idx]+str(col_idx+1)
    expected={s['sample_name'] for s in samples}
    missing=expected-set(mapping)
    if missing:return {}, f'Miniprep worksheet missing {len(missing)} colonies: {", ".join(sorted(missing)[:6])}. Fix the layout worksheet.'
    return mapping,None


def write_html_report(outpath, chain, plate_id, data_name, concise, samples, rows, miniprep_mapping=None, miniprep_warning=None):
    """One printable offline report: action matrix plus sequencing plate locator."""
    from html import escape as h
    from collections import Counter
    groups = ['USE','USE - WRONG WELL','PURIFY','PURIFY - WRONG WELL',
              'SCREEN MORE COLONIES','RETRANSFORM GIBSON','CHECK FILES']
    labels = {
       'USE': 'Ready to use', 'USE - WRONG WELL':'Correct insert found in an unexpected colony',
       'PURIFY':'Purify and screen', 'PURIFY - WRONG WELL':'Correct insert found in an unexpected colony — purify',
       'SCREEN MORE COLONIES':'Screen more colonies',
       'RETRANSFORM GIBSON':'Retransform Gibson', 'CHECK FILES':'Check sequencing inputs'
    }
    instructions = {
       'USE':'Pick the indicated colony from its patch plate, confirm patch identity, and expand for maxiprep.',
       'USE - WRONG WELL':'This colony contains the correct insert but was assigned to a different construct. Retrieve it using its original colony identifier on the patch plate, then confirm its identity by sequencing an independently recovered isolate before expansion.',
       'PURIFY':'Re-isolate a single clone from the indicated sample and verify by sequencing.',
       'PURIFY - WRONG WELL':'The correct insert was found in a colony assigned to another construct, but its sequence purity is insufficient. Locate it by its original colony identifier, isolate a single clone, and confirm by sequencing before proceeding.',
       'SCREEN MORE COLONIES':'Screen additional colonies from this transformation; if none remain, retransform the Gibson assembly.',
       'RETRANSFORM GIBSON':'Retransform the existing Gibson assembly and pick new colonies.',
       'CHECK FILES':'Correct missing, mismatched or unreadable sequencing files, then rerun.'
    }
    cssclasses = {'USE':'use','USE - WRONG WELL':'wrong','PURIFY':'purify',
                  'PURIFY - WRONG WELL':'wrong','SCREEN MORE COLONIES':'screen',
                  'RETRANSFORM GIBSON':'redo','CHECK FILES':'check'}
    counts=Counter(r['ACTION'] for r in concise)
    use_n=sum(counts[x] for x in ('USE','USE - WRONG WELL'))
    pur_n=sum(counts[x] for x in ('PURIFY','PURIFY - WRONG WELL'))
    work_n=len(concise)-use_n-pur_n
    src_by_num={r['sample_number']:r for r in rows}
    selected={int(r['SEQ_SAMPLE_NUMBER']):r for r in concise if str(r['SEQ_SAMPLE_NUMBER']).strip()}
    def fmt(v):
        if v is None or str(v).lower() in ('nan','none'):return '—'
        return h(str(v)) if str(v) else '—'
    def idt(r):return f'<strong>{fmt(r["IDT_WELL"])}</strong><small>{fmt(r["INSERT"])}</small>'
    action_sections=[]
    for act in groups:
        group=[r for r in concise if r['ACTION']==act]
        if not group:continue
        lines=[]
        for r in group:
            chosen=r['PICK_COLONY']
            detail=[]
            if chosen:
                detail.append(f'<b>{fmt(chosen)}</b>')
                if r['SEQ_SAMPLE_NUMBER']:
                    idx=int(r['SEQ_SAMPLE_NUMBER'])-1
                    detail.append(f'sequencing plate {"ABCDEFGH"[idx//12]}{idx%12+1}')
                if r['MEDIAN_CONSENSUS_%']!='':detail.append(f'{fmt(r["MEDIAN_CONSENSUS_%"])}% insert consensus')
                if r['SYNONYMOUS_MUTATIONS']: detail.append(f'silent changes: {fmt(r["SYNONYMOUS_MUTATIONS"])}')
            else:
                detail.append('No accepted colony')
            alternatives=[]
            if r['ALL_GOOD_COLONIES']:
                alternatives.append(f'Usable: {fmt(r["ALL_GOOD_COLONIES"])}')
            if r['ALL_PURIFY_COLONIES']:
                alternatives.append(f'Purification candidates: {fmt(r["ALL_PURIFY_COLONIES"])}')
            caveat=''
            if chosen and r['PLASMID_INTEGRITY'] != 'intact':
                info = r['INTEGRITY_FLAGGED_WINDOWS'] or 'QC not assessed'
                caveat=f'<div class="warning">Backbone QC follow-up before proceeding: {fmt(info)}</div>'
            if 'WRONG WELL' in act:
                origin=r['ACTUAL_ORIGIN_WELL']
                caveat += f'<div class="warning">Recorded original source: {fmt(origin)}. Verify patch provenance.</div>'
            lines.append('<tr><td>'+idt(r)+'</td><td>'+ ' · '.join(detail)+ (f'<div class="minor">{" | ".join(alternatives)}</div>' if alternatives else '')+caveat+'</td></tr>')
        action_sections.append(f'<section class="section"><h2><span class="dot {cssclasses[act]}"></span>{labels[act]} <span class="n">{len(group)}</span></h2><p class="instruction">{instructions[act]}</p><table><thead><tr><th>eBlock / original IDT well</th><th>Colony and decision details</th></tr></thead><tbody>{"".join(lines)}</tbody></table></section>')
    def well_html(position, name, category, subtitle='', tooltip=''):
        return (f'<div class="well {category}" title="{h(tooltip or name)}">'
                f'<span class="well-pos">{h(position)}</span><strong>{h(name)}</strong>'
                f'<small>{h(subtitle)}</small></div>')

    def outcome(r):
        if not r or r.get('decision')=='review':return 'check','Check data'
        decision=r.get('decision','')
        if decision=='failed': return 'failed','Seq failed'
        if decision in ('wrong_well_ready','wrong_well_purify'):return 'wrong','Other insert'
        if decision=='ready':return 'use','Correct, pure'
        if decision=='purify':return 'purify','Purify'
        return 'bad','Too mixed / wrong'

    sequencing_map={}
    for source in samples:
        pos='ABCDEFGH'[(source['sample_number']-1)//12]+str((source['sample_number']-1)%12+1)
        sequencing_map[pos]=source['sample_name']

    def render_colony_plate(mapping):
        output=[]
        for i in range(96):
            pos='ABCDEFGH'[i//12]+str(i%12+1)
            colony=mapping.get(pos)
            if colony is None:
                output.append(well_html(pos,'','empty'))
                continue
            r=by_colony.get(colony)
            cls,sub=outcome(r)
            observed=r.get('best_candidate_well','') if r else ''
            if cls=='wrong' and observed: sub='Insert '+str(observed)
            title=f'{pos}: {colony}; {sub}; {r.get("notes", "") if r else "not sequenced"}'
            output.append(well_html(pos,colony,cls,sub,title))
        return ''.join(output)

    by_colony={r['sample_name']:r for r in rows}
    sequencing_plate=render_colony_plate(sequencing_map)
    miniprep_plate=render_colony_plate({v:k for k,v in (miniprep_mapping or {}).items()}) if not miniprep_warning else ''
    assembly_map={r['ASSEMBLY_WELL']:r for r in concise}
    gibson_cells=[]
    for i in range(96):
        pos='ABCDEFGH'[i//12]+str(i%12+1)
        r=assembly_map.get(pos)
        if not r:
            gibson_cells.append(well_html(pos,'','empty'))
            continue
        action=r['ACTION'];cls=cssclasses.get(action,'check')
        pick=r.get('PICK_COLONY') or ''
        short={'USE':'Use','USE - WRONG WELL':'Use (other colony)',
               'PURIFY':'Purify','PURIFY - WRONG WELL':'Purify (other)',
               'SCREEN MORE COLONIES':'Pick more','RETRANSFORM GIBSON':'Retransform',
               'CHECK FILES':'Check files'}.get(action,action)
        gibson_cells.append(well_html(pos,r['IDT_WELL'],cls,pick or short,
                                      f'Assembly {pos}; eBlock {r["IDT_WELL"]}: {r["INSERT"]}; {action}; preferred colony {pick}'))
    gibson_plate=''.join(gibson_cells)
    headings=''.join(f'<span>{i}</span>' for i in range(1,13))
    miniprep_content=(f'<div class="plate-grid">{miniprep_plate}</div>' if not miniprep_warning
                      else '<p class="warning">'+h(miniprep_warning)+'</p>')
    gallery=f'''<section class="plate-section"><h2>Plate maps</h2>
    <p class="instruction">Switch between the physical plate layouts. Labels give the original colony identifier (or IDT well on the Gibson plate); hover over a well for more information.</p>
    <div class="tabs" role="tablist" aria-label="Plate maps"><button type="button" class="active" onclick="showPlate('seq',this)">Sequencing</button><button type="button" onclick="showPlate('mini',this)">Miniprep</button><button type="button" onclick="showPlate('gib',this)">Gibson assembly</button></div>
    <div id="plate-seq" class="plate-panel"><p class="instruction">Sequencing sample positions (packed by sample number).</p><div class="plate-scroll"><div class="plate-heading">{headings}</div><div class="plate-grid">{sequencing_plate}</div></div></div>
    <div id="plate-mini" class="plate-panel" hidden><p class="instruction">Original miniprep positions from the workbook. These are NOT sequencing positions.</p><div class="plate-scroll"><div class="plate-heading">{headings}</div>{miniprep_content}</div></div>
    <div id="plate-gib" class="plate-panel" hidden><p class="instruction">Assembly positions after the documented eBlock well relocations. Each well shows its original IDT well and best available colony or next action.</p><div class="plate-scroll"><div class="plate-heading">{headings}</div><div class="plate-grid">{gibson_plate}</div></div></div>
    <div class="legend"><span><i class="dot use"></i>Correct, pure / ready</span><span><i class="dot wrong"></i>Correct insert in unexpected colony</span><span><i class="dot purify"></i>Correct, purify</span><span><i class="dot bad"></i>Mixed / incorrect</span><span><i class="dot failed"></i>Sequencing failed</span><span><i class="dot screen"></i>Screen more</span><span><i class="dot redo"></i>Retransform</span><span><i class="dot check"></i>Check input</span></div></section>'''
    stat=f'''<div class="stats"><div><strong>{len(concise)}</strong><span>Intended inserts</span></div><div><strong>{use_n}</strong><span>Ready-to-use candidates</span></div><div><strong>{pur_n}</strong><span>Need purification</span></div><div><strong>{work_n}</strong><span>Need additional work</span></div></div>'''
    html=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{h(chain)} sequencing — action report</title>
<style>
:root{{--ink:#1e293b;--muted:#617084;--line:#dce3eb;--bg:#f6f8fb}}
*{{box-sizing:border-box}}body{{font-family:system-ui,-apple-system,'Segoe UI',sans-serif;color:var(--ink);background:var(--bg);margin:0;line-height:1.42}}main{{max-width:1180px;margin:auto;padding:30px 24px 80px}}header{{margin:0 0 24px}}h1{{font-size:30px;letter-spacing:-.03em;margin:0}}.meta{{color:var(--muted);margin:6px 0 0}}.stats{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:23px 0 30px}}.stats>div{{background:white;border:1px solid var(--line);border-radius:10px;padding:16px}}.stats strong{{display:block;font-size:30px}}.stats span{{font-size:13px;color:var(--muted)}}.section,.plate-section{{background:white;border:1px solid var(--line);border-radius:12px;margin:18px 0;padding:20px 22px;break-inside:avoid}}h2{{font-size:19px;margin:0;display:flex;align-items:center;gap:10px}}.n{{font-size:13px;font-weight:600;background:#f0f3f7;border-radius:50px;padding:2px 9px}}.instruction{{margin:8px 0 15px;color:#4c5867;font-size:14px}}table{{border-collapse:collapse;width:100%;table-layout:fixed}}td,th{{text-align:left;padding:11px 10px;border-top:1px solid var(--line);vertical-align:top;font-size:14px;overflow-wrap:anywhere}}th{{color:var(--muted);font-weight:600;font-size:12px;text-transform:uppercase}}th:first-child,td:first-child{{width:35%}}td small{{display:block;color:var(--muted);font-size:12px;margin-top:2px}}.minor{{color:#526276;font-size:12px;margin-top:4px}}.warning{{color:#8c3c10;background:#fff3e8;border-left:3px solid #e8a04c;padding:6px 8px;margin-top:6px;font-size:12px}}.dot{{display:inline-block;width:11px;height:11px;border-radius:3px;flex-shrink:0}}.use{{background:#c7e8d3}}.wrong{{background:#ffd8a9}}.purify{{background:#f5e3a8}}.screen{{background:#d7e7fc}}.redo{{background:#e9ced3}}.check{{background:#eadff9}}.other{{background:#eff1f4}}.bad{{background:#d3d8e0}}.failed{{background:#697383;color:white}}.empty{{background:#f8fafc}}.plate-scroll{{overflow-x:auto;padding-bottom:5px}}.plate-grid{{display:grid;grid-template-columns:repeat(12,minmax(66px,1fr));gap:4px;min-width:830px}}.plate-heading{{display:grid;grid-template-columns:repeat(12,1fr);min-width:830px;gap:4px;text-align:center;margin-bottom:5px;color:var(--muted);font-size:12px;font-weight:600}}.well{{min-height:72px;border-radius:6px;border:1px solid #cdd5dd;padding:6px;text-align:center;position:relative}}.well-pos{{font-size:10px;display:block;color:#526276;text-align:left}}.well strong{{font-size:11px;display:block;margin-top:4px}}.well small{{font-size:10px;display:block;color:#526276}}.legend{{display:flex;flex-wrap:wrap;gap:14px;font-size:12px;margin:14px 0 0;color:#45556a}}.legend span{{display:flex;align-items:center;gap:5px}}.tabs{{display:flex;gap:7px;margin:15px 0 20px}}.tabs button{{font:inherit;background:#e8eef6;border:1px solid #ccd6e1;border-radius:8px;padding:9px 15px;cursor:pointer}}.tabs button.active{{background:#245783;color:white}}.plate-panel[hidden]{{display:none}}.footer{{font-size:12px;color:var(--muted);margin-top:20px}}@media(max-width:700px){{main{{padding:16px}}.stats{{grid-template-columns:repeat(2,1fr)}}.section{{padding:14px}}th:first-child,td:first-child{{width:40%}}}}@media print{{body{{background:white}}main{{max-width:none;padding:0}}.section,.plate-section,.stats>div{{border-color:#ccc}}.plate-section{{break-before:page}}.well,.dot{{print-color-adjust:exact;-webkit-print-color-adjust:exact}}}}
</style></head><body><main><header><h1>{h(chain)} plasmid sequencing</h1><p class="meta">Plate {h(plate_id)} · {h(data_name)} · {len(samples)} sequencing samples · One decision per eBlock</p></header>{stat}
<div class="section"><h2>Action plan</h2><p class="instruction">Use the preferred colony for each eBlock. An insert with a usable colony does not also require purification. Orange source discrepancies require identity confirmation. Backbone warnings override a ready-to-use designation until resolved.</p></div>
{''.join(action_sections)}
{gallery}
<p class="footer">Sequence acceptance: full-length expected eBlock, approved synonymous changes and TCRb junction substitution; median consensus &gt;{PURE_MEDIAN_THRESHOLD}% for use or &gt;{USABLE_MEDIAN_THRESHOLD}% for purification. Backbone homogeneity screen: {BACKBONE_WINDOW_BP}-bp windows outside insert, median &ge;{BACKBONE_MIN_MEDIAN_CONSENSUS}%. “Intact” denotes a passing homogeneity screen, not independent structural validation. Provider Fail is treated as no plasmid.</p></main><script>function showPlate(id,button){{document.querySelectorAll('.plate-panel').forEach(x=>x.hidden=true);document.getElementById('plate-'+id).hidden=false;document.querySelectorAll('.tabs button').forEach(b=>b.classList.remove('active'));button.classList.add('active');}}</script></body></html>'''
    outpath.write_text(html, encoding='utf-8')

def write_actionable_reports(data_dir, all_eblocks, layouts, chain_map, diagnostics=False):
    plate_id,chain=select_chain(data_dir,chain_map)
    current={w:(name,ref) for w,(name,ref) in all_eblocks.items()
             if ('TCRa' if w[0] in 'ABCD' else 'TCRb')==chain}
    samples=sorted((x for x in layouts if x['plate']==chain),key=lambda x:x['sample_number'])
    outdir=data_dir/'results';outdir.mkdir(exist_ok=True)
    detaildir=outdir/'diagnostics'
    if diagnostics: detaildir.mkdir(exist_ok=True)
    summary=load_sequencing_summary(data_dir,plate_id)
    fasta_dir=locate_folder(data_dir,'FASTA_FILES')
    per_dir=locate_folder(data_dir,'PER_BASE_DATA')
    aln=make_aligner()
    rows=[]
    candidate_to_well=defaultdict(list)
    samples_by_expected=defaultdict(list)
    print(f'  Analyzing {len(samples)} sequencing samples for {chain} ({plate_id})...', flush=True)
    for progress_index, sm in enumerate(samples, 1):
        sn=sm['sample_name'];n=sm['sample_number'];expected=sm['idt_well']
        rec=summary.get(n) if summary else None
        result=dict(sample_number=n,sample_name=sn,original_assembly_well=sm['assembly_well'],
                    original_idt_well=expected, sequencing_status=(rec or {}).get('status','unknown'),
                    decision='review',notes='')
        if rec and rec['name']!=sn:
            result.update(decision='review',notes='Summary sample identifier disagrees with layout')
            rows.append(result);samples_by_expected[expected].append(result)
            print(f'    [{progress_index}/{len(samples)}] {sn}: {result["decision"]}', flush=True)
            continue
        if rec and rec['status'].lower()=='fail':
            result.update(decision='failed',notes='Provider Fail interpreted as plasmid absent')
            rows.append(result);samples_by_expected[expected].append(result)
            print(f'    [{progress_index}/{len(samples)}] {sn}: {result["decision"]}', flush=True)
            continue
        f=find_sample_file(fasta_dir,n,sn,'.fasta') or find_sample_file(fasta_dir,n,sn,'.fa')
        b=find_sample_file(per_dir,n,sn,'.csv')
        result.update(fasta_file=str(f.relative_to(data_dir)) if f else '',
                      per_base_file=str(b.relative_to(data_dir)) if b else '')
        if not f or not b:
            result.update(decision='review',notes='Missing FASTA or per-base file')
            rows.append(result);samples_by_expected[expected].append(result)
            print(f'    [{progress_index}/{len(samples)}] {sn}: {result["decision"]}', flush=True)
            continue
        try:
            seq=get_seq(f)
            consensus,positions,percent=get_bases(b)
            exact=_exact_refs(seq,all_eblocks,aln)
            if expected not in exact and matching_reference(expected,all_eblocks[expected][1],seq,aln):exact.append(expected)
            same_chain=[w for w in exact if w in current]
            other_chain=[w for w in exact if w not in current]
            result['exact_reference_wells']='; '.join(exact)
            result['exact_other_chain_wells']='; '.join(other_chain)
            if same_chain:
                # Determine purity using one full-length EXACT reference sequence.
                chosen=expected if expected in same_chain else same_chain[0]
                original_reference=all_eblocks[chosen][1]
                approved=classify_reference_match(chosen,original_reference,seq,aln)
                matched_reference=approved['sequence']
                uses_variant=approved['allowed_junction_variant']
                median,coverage=_consensus_for_reference(aln,matched_reference,consensus,percent)
                integrity, affected, windows=assess_backbone_integrity(
                    aln,matched_reference,consensus,positions,percent)
                result.update(plasmid_integrity=integrity,
                              integrity_flagged_windows=affected,
                              backbone_windows_evaluated=len(windows))
                result.update(allowed_junction_variant=uses_variant,
                              synonymous_mismatches=approved['synonymous_mismatches'],
                              synonymous_changes=approved['synonymous_changes'],
                              best_candidate_well=chosen,
                              best_candidate_name=all_eblocks[chosen][0],
                              best_candidate_identity_pct=100.0,
                              best_candidate_coverage_pct=100.0,
                              best_candidate_mismatches=0,best_candidate_insertions_bp=0,
                              best_candidate_deletions_bp=0,median_consensus_pct=median,
                              per_base_coverage_pct=coverage)
                if median is None:
                    result.update(decision='review',notes='Exact insert in FASTA, but per-base insert location incomplete')
                elif median>PURE_MEDIAN_THRESHOLD:
                    result.update(decision='ready' if expected in same_chain else 'wrong_well_ready',
                                  notes='Full-length accepted insert; synonymous: '+approved['synonymous_changes']+('; allowed TCRb junction G->A' if uses_variant else ''))
                elif median>USABLE_MEDIAN_THRESHOLD:
                    result.update(decision='purify' if expected in same_chain else 'wrong_well_purify',
                                  notes='Full-length accepted insert; synonymous: '+approved['synonymous_changes']+('; allowed TCRb junction G->A' if uses_variant else ''))
                else:
                    result.update(decision='inadequate',notes='Exact insert but consensus too mixed')
                if result['decision'] in ('ready','wrong_well_ready','purify','wrong_well_purify'):
                    for well in same_chain:
                        # Identical eBlocks cannot establish a swapped colony.
                        # Attribute such a sample to its intended well only.
                        if well == expected or (expected not in same_chain and all_eblocks[well][1] != all_eblocks[expected][1]):
                            candidate_to_well[well].append(result)
            else:
                # Search approximate candidates for review only; never grant a 'good' call.
                candidates=_rank_kmers(seq,all_eblocks)
                if expected not in candidates:candidates.append(expected)
                scored=[]
                for well in candidates:
                    stats=evaluate(aln,all_eblocks[well][1],seq)
                    scored.append((stats['coverage_pct']>=MIN_EXPECTED_COVERAGE-1e-7,
                                   stats['identity_pct'],stats['coverage_pct'],well,stats))
                scored.sort(key=lambda t:(t[0],t[1],t[2]),reverse=True)
                chosen=scored[0]
                match=chosen[4]
                median,coverage=_consensus_for_reference(aln,all_eblocks[chosen[3]][1],consensus,percent)
                result.update(best_candidate_well=chosen[3],
                              best_candidate_name=all_eblocks[chosen[3]][0],
                              best_candidate_identity_pct=round(match['identity_pct'],4),
                              best_candidate_coverage_pct=round(match['coverage_pct'],4),
                              best_candidate_mismatches=match['mismatches'],
                              best_candidate_insertions_bp=match['insertions_bp'],
                              best_candidate_deletions_bp=match['deletions_bp'],
                              median_consensus_pct=median,per_base_coverage_pct=coverage,
                              decision='needs_review' if match['identity_pct']>=DIAGNOSTIC_MIN_IDENTITY
                                       else 'inadequate',
                              notes='Non-exact match; cannot qualify for maxiprep or purification; inspect diagnostic candidates')
                # keep complete ranked evidence quiet by default
        except Exception as ex:
            result.update(decision='review',notes=f'Analysis error: {type(ex).__name__}: {ex}')
        rows.append(result);samples_by_expected[expected].append(result)
        print(f'    [{progress_index}/{len(samples)}] {sn}: {result["decision"]}', flush=True)
    if diagnostics: output_csv(detaildir/'colony_diagnostics.csv',DETAIL_FIELDS,rows)
    action_rows=[]
    for well,(name,ref) in current.items():
        expected_samples=samples_by_expected.get(well,[])
        potential=candidate_to_well.get(well,[])
        # Exclude same-name duplicates from being flagged wrong-well. They cannot
        # be distinguished by the insert; retain provenance caveat.
        def rank(r):
            return (r['decision'] in ('ready','wrong_well_ready'),
                    r.get('plasmid_integrity')=='intact',
                    r['original_idt_well']==well,
                    float(r.get('median_consensus_pct') or 0))
        potential.sort(key=rank,reverse=True)
        best=potential[0] if potential else None
        pure_candidates=[r for r in potential if r['decision'] in ('ready','wrong_well_ready')]
        purify_candidates=[r for r in potential if r['decision'] in ('purify','wrong_well_purify')]
        if best:
            correct_origin=best['original_idt_well']==well
            is_pure=best['decision'] in ('ready','wrong_well_ready')
            if is_pure:
                action='USE' if correct_origin else 'USE - WRONG WELL'
            else:
                action='PURIFY' if correct_origin else 'PURIFY - WRONG WELL'
            reason=best.get('notes','')
            if best.get('plasmid_integrity')=='follow up':
                reason+='; BACKBONE FOLLOW UP BEFORE USE: '+best.get('integrity_flagged_windows','')
            elif best.get('plasmid_integrity')!='intact':
                reason+='; BACKBONE NOT ASSESSED: inspect per-base data'
            # Identical reference constructs are not distinguishable by insert sequence.
            same_seq=[w for w,(nm,s) in all_eblocks.items() if s==ref]
            if len(same_seq)>1:reason+='; same-sequence references: '+', '.join(same_seq)
        elif not expected_samples:
            action='SCREEN MORE COLONIES';reason='NO SEQUENCING SAMPLE: screen saved colonies; if none exist, retransform Gibson assembly'
        elif any(r['decision']=='review' for r in expected_samples):
            action='CHECK FILES';reason='Sequencing input or sample mapping missing/inconsistent; correct files or map then rerun'
        elif any(r['decision']=='needs_review' for r in expected_samples):
            action='SCREEN MORE COLONIES'
            reason='No acceptable full-length sequence in this set; screen other colonies from this transformation; if exhausted, retransform Gibson assembly'
        else:
            action='RETRANSFORM GIBSON'
            reason='No acceptable candidate; failed sequencing counts as plasmid absent. Repeat assembly transformation and screen new colonies'
        action_rows.append(dict(chain=chain,idt_well=well,
                 assembly_well=next((k for k,v in RELOCATIONS.items() if v==well),well),
                 insert_name=name,action=action,selected_sample=best['sample_name'] if best else '',
                 sequencing_well=best['sample_number'] if best else '',
                 origin_assembly_well=best['original_assembly_well'] if best else '',
                 median_consensus_pct=best.get('median_consensus_pct','') if best else '',
                 identity_pct=best.get('best_candidate_identity_pct','') if best else '',
                 coverage_pct=best.get('best_candidate_coverage_pct','') if best else '',
                 good_colonies='; '.join(r['sample_name'] for r in pure_candidates),
                 purify_colonies='; '.join(r['sample_name'] for r in purify_candidates),
                 other_candidate_samples='; '.join(r['sample_name'] for r in potential[1:]),
                 synonymous_mismatches=best.get('synonymous_mismatches',0) if best else '',
                 synonymous_changes=best.get('synonymous_changes','') if best else '',
                 plasmid_integrity=best.get('plasmid_integrity','not assessed') if best else 'not assessed',
                 integrity_flagged_windows=best.get('integrity_flagged_windows','') if best else '',
                 reason=reason))
    # One compact, filtered CSV per sequencing plate. One row per intended eBlock.
    from collections import Counter
    concise=[]
    for r in action_rows:
        concise.append({
            'ACTION':r['action'], 'IDT_WELL':r['idt_well'],
            'INSERT':r['insert_name'], 'ASSEMBLY_WELL':r['assembly_well'],
            'PICK_COLONY':r['selected_sample'],
            'SEQ_SAMPLE_NUMBER':r['sequencing_well'],
            'ACTUAL_ORIGIN_WELL':r['origin_assembly_well'],
            'MEDIAN_CONSENSUS_%':r['median_consensus_pct'],
            'SYNONYMOUS_MUTATIONS':r['synonymous_changes'],
            'PLASMID_INTEGRITY':r['plasmid_integrity'],
            'INTEGRITY_FLAGGED_WINDOWS':r['integrity_flagged_windows'],
            'ALL_GOOD_COLONIES':r['good_colonies'],
            'ALL_PURIFY_COLONIES':r['purify_colonies'],
            'WHAT_TO_DO':r['reason'],
        })
    fields=list(concise[0]) if concise else ['ACTION','IDT_WELL','INSERT']
    order=['USE','USE - WRONG WELL','PURIFY','PURIFY - WRONG WELL',
           'SCREEN MORE COLONIES','RETRANSFORM GIBSON','CHECK FILES']
    concise.sort(key=lambda r:(order.index(r['ACTION']) if r['ACTION'] in order else 99,r['IDT_WELL']))
    filename=outdir/f'{chain}_ACTIONS.csv'
    output_csv(filename,fields,concise)
    report_path=outdir/f'{chain}_REPORT.html'
    miniprep_mapping,miniprep_warning=read_miniprep_mapping(Path(__file__).resolve().parent / "layouts",chain,samples)
    write_html_report(report_path,chain,plate_id,data_dir.name,concise,samples,rows,miniprep_mapping,miniprep_warning)
    totals=Counter(r['ACTION'] for r in concise)
    # Remove only the obsolete output files created by previous versions.
    oldfiles=['colony_results.csv','insert_summary.csv','verified_good_and_pure.csv',
              'action_required.csv','summary.txt','final_decisions_by_insert.csv',
              'good_to_use.csv','retransform_to_purify.csv',
              'retransform_from_gibson.csv','correct_sequence_wrong_well.csv',
              'screen_more_or_review.csv','FINAL_README.txt',
              f'{chain}_WET_LAB_ACTIONS.csv', f'{chain}_SUMMARY.txt']
    for name in oldfiles:
        (outdir/name).unlink(missing_ok=True)
    if not diagnostics and detaildir.is_dir():
        import shutil
        shutil.rmtree(detaildir)
    print(f'  Written: {filename}')
    print(f'  Written: {report_path}')
    for status in order:
        if totals[status]: print(f'    {status}: {totals[status]} eBlocks requiring this action')
    n_purify=sum(bool(r['purify_colonies']) for r in action_rows)
    n_purify_only=sum(r['action'].startswith('PURIFY') for r in action_rows)
    print(f'    Note: {n_purify} eBlocks have >=1 purification candidate, '
          f'but only {n_purify_only} REQUIRE purification because a good clone '
          'takes precedence where one exists.')

if __name__=='__main__':main()
