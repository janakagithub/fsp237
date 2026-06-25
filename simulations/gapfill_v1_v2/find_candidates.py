#!/usr/bin/env python3
"""Phase 2 + Phase 3 -- find Colletotrichum gene candidates for each gapfill
reaction in Version 1, then BLAST against the FSP237 proteome.

Inputs
  - simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version1_noGenes.json
  - gpr-update/C_higgensium.gbff               (C. higginsianum IMI 349063 RefSeq)
  - gpr-update/blast_db/fsp237                 (FSP237 blastp DB; existing)
  - gpr-update/Csublineola_reference_plus_novel_classu.proteins.fa

Outputs (under simulations/gapfill_v1_v2/)
  - candidates/v1_query_candidates.faa         (FASTA: chosen C. higg proteins)
  - candidates/v1_query_provenance.tsv         (mapping: rxn_id, locus_tag, product, why-this-protein)
  - blast/v1_blast_results.tsv                 (per-query top FSP237 hits)
  - blast/v1_rxn_to_gene_mapping.tsv           (final picks per rxn_id with metrics + confidence)

Logic per gap-fill reaction:
  1. Score every C. higginsianum protein against the reaction's keyword set
     (product keywords + EC numbers parsed from /product= and /EC_number=).
     Keep the top 1-2 candidates.
  2. Pull their protein sequences from C_higgensium.gbff.
  3. BLAST each candidate against fsp237 blast_db.
  4. Apply quality thresholds (pident >= 30, qcov >= 50, e-value <= 1e-10);
     classify each rxn as confident / weak / no_hit and emit a recommended
     gene_* assignment for Phase 4.
"""
import os
import re
import subprocess
import sys
import textwrap

# ---- spec table: reactions vs keywords / EC / source-strain hint --------
# Tied to the V1 gapfill list. Keywords are matched against /product= text;
# EC numbers are exact-matched against /EC_number= when present.
GAPFILL_TARGETS = [
    # (rxn_id, label, ec_set, name_keywords, source_strain_hint)

    # --- A1. Acyl-CoA dehydrogenases (peroxisomal, chain-length-specific) ---
    # In fungi a single Fox1 / acyl-CoA oxidase (POX1) handles a broad chain
    # range; the same CH63R protein is therefore the right candidate for all
    # three chain lengths (C4/C6/C8), but we score per-rxn so each call lands
    # the right hit if the genome encodes specific isoforms.
    ('rxn02679_x0', 'Octanoyl-CoA:O2 2-oxidoreductase (C8)',
        {'1.3.3.6','1.3.99.3'},
        ['acyl-coa oxidase', 'acyl-coenzyme a oxidase', 'pox1', 'pox', 'fox1'],
        'C. higginsianum'),
    ('rxn03251_x0', 'Hexanoyl-CoA 2,3-oxidoreductase (C6)',
        {'1.3.3.6','1.3.99.3'},
        ['acyl-coa oxidase', 'acyl-coenzyme a oxidase', 'pox1', 'pox', 'fox1'],
        'C. higginsianum'),
    ('rxn00868_x0', 'Butanoyl-CoA 2,3-oxidoreductase (C4)',
        {'1.3.8.1','1.3.99.2'},
        ['butyryl-coa dehydrogenase', 'short chain', 'acyl-coa dehydrogenase'],
        'C. higginsianum'),

    # --- A2. Enoyl-CoA hydratases (fungal Fox2 / MFP handles all chains) ---
    ('rxn03247_x0', '(S)-Hydroxyoctanoyl-CoA hydratase (C8)',
        {'4.2.1.17','4.2.1.74'},
        ['enoyl-coa hydratase', 'hydroxyacyl-coa', 'multifunctional',
         'multifunctional enzyme', 'mfp', 'fox2', 'isomerase'],
        'C. higginsianum'),
    ('rxn03250_x0', '(S)-Hydroxyhexanoyl-CoA hydratase (C6)',
        {'4.2.1.17','4.2.1.74'},
        ['enoyl-coa hydratase', 'hydroxyacyl-coa', 'multifunctional',
         'multifunctional enzyme', 'mfp', 'fox2'],
        'C. higginsianum'),
    ('rxn02167_x0', '(S)-3-Hydroxybutanoyl-CoA hydratase (C4)',
        {'4.2.1.55','4.2.1.17'},
        ['enoyl-coa hydratase', 'crotonase', 'crotonyl-coa', 'hydroxybutyryl-coa'],
        'C. higginsianum'),

    # --- A3. 3-OH-acyl-CoA dehydrogenases (fungal MFP/Fox2 - same enzyme) ---
    ('rxn03246_x0', '(S)-hydroxyoctanoyl-CoA dehydrogenase (C8)',
        {'1.1.1.35'},
        ['3-hydroxyacyl-coa dehydrogenase', 'multifunctional', 'fox2', 'mfp'],
        'C. higginsianum'),
    ('rxn03249_x0', '(S)-hydroxyhexanoyl-CoA dehydrogenase (C6)',
        {'1.1.1.35'},
        ['3-hydroxyacyl-coa dehydrogenase', 'multifunctional', 'fox2', 'mfp'],
        'C. higginsianum'),
    ('rxn03861_x0', '(S)-3-hydroxybutanoyl-CoA / acetoacetyl-CoA reductase (C4)',
        {'1.1.1.36','1.1.1.35'},
        ['acetoacetyl-coa reductase', '3-hydroxybutyryl-coa dehydrogenase',
         'aldo-keto reductase'],
        'C. higginsianum'),

    # --- A4. Chain-shortening 3-ketoacyl-CoA thiolases (Pot1 family) ---
    # In yeast/fungi 1 or 2 thiolases handle most chain lengths.
    ('rxn02804_x0', '3-ketoacyl-CoA thiolase (C16->C14)',
        {'2.3.1.16','2.3.1.9'},
        ['ketoacyl-coa thiolase', '3-ketoacyl-coa thiolase', 'acyltransferase',
         'pot1', 'fox3'],
        'C. higginsianum'),
    ('rxn06510_x0', '3-ketoacyl-CoA thiolase (C14->C12)',
        {'2.3.1.16'},
        ['ketoacyl-coa thiolase', 'acyltransferase', 'pot1', 'fox3'],
        'C. higginsianum'),
    ('rxn03243_x0', '3-ketoacyl-CoA thiolase (C12->C10)',
        {'2.3.1.16'},
        ['ketoacyl-coa thiolase', 'acyltransferase', 'pot1', 'fox3'],
        'C. higginsianum'),
    ('rxn02680_x0', '3-ketoacyl-CoA thiolase (C10->C8)',
        {'2.3.1.16'},
        ['ketoacyl-coa thiolase', 'acyltransferase', 'pot1', 'fox3'],
        'C. higginsianum'),
    ('rxn03248_x0', '3-ketoacyl-CoA thiolase (C8->C6)',
        {'2.3.1.16'},
        ['ketoacyl-coa thiolase', 'acyltransferase', 'pot1', 'fox3'],
        'C. higginsianum'),
    ('rxn00874_x0', '3-ketoacyl-CoA thiolase (C6->C4)',
        {'2.3.1.16'},
        ['ketoacyl-coa thiolase', 'acyltransferase', 'pot1', 'fox3'],
        'C. higginsianum'),
    ('rxn00178_x0', 'Acetyl-CoA C-acetyltransferase (C4->2 AcCoA)',
        {'2.3.1.9'},
        ['acetyl-coa c-acetyltransferase', 'acetoacetyl-coa thiolase',
         'acetyl-coa acetyltransferase', 'erg10', 'thiolase'],
        'C. higginsianum'),

    # --- A5. Peroxisomal cofactor shuttles (NEW transport stubs) ---
    # These don't necessarily have a single dedicated transporter gene; the
    # Ant1 ortholog covers ATP/AMP, and yeast PXA1/2 covers fatty acyl-CoA.
    # NAD/CoA/PPi/H+ likely move via porin-like or undefined routes; assign
    # by similarity to Ant1 / generic peroxisomal-membrane MFS.
    ('tx_atp_xc', 'ATP/AMP antiport (Ant1-like, peroxisomal)',
        set(),
        ['ant1', 'adenine nucleotide carrier', 'peroxisomal adenine nucleotide'],
        'C. higginsianum'),
    ('tx_ppi_xc', 'PPi peroxisomal export',
        set(),
        ['pyrophosphatase', 'inorganic pyrophosphatase', 'pppi'],
        'C. higginsianum'),
    ('tx_coa_xc', 'CoA peroxisomal equilibration',
        set(),
        ['coa transporter', 'peroxisomal coa', 'mfs'],
        'C. higginsianum'),
    ('tx_nad_xc', 'NAD peroxisomal equilibration',
        set(),
        ['nadhx', 'pyridine nucleotide', 'nad transporter', 'mitochondrial carrier'],
        'C. higginsianum'),
    ('tx_nadh_xc', 'NADH peroxisomal equilibration',
        set(),
        ['nadhx', 'pyridine nucleotide', 'nadh transporter', 'mitochondrial carrier'],
        'C. higginsianum'),
    ('tx_h_xc', 'H+ peroxisomal equilibration',
        set(),
        ['proton transporter', 'peroxisome', 'h+'],
        'C. higginsianum'),
    # ---- L-arabinose pathway ----
    ('rxn01391_c0', 'L-arabitol:NAD+ 4-oxidoreductase (L-xylulose-forming)',
        {'1.1.1.12'}, ['l-arabitol', 'arabinitol 4-dehydrogenase', 'arabinitol',
                        'lad1', 'l-arabitol dehydrogenase'], 'C. higginsianum'),
    ('rxn33066_c0', 'L-xylulose reductase',
        {'1.1.1.10'}, ['l-xylulose reductase', 'xylulose reductase', 'lxr',
                        'lxr1', 'lxr4'], 'C. higginsianum'),
    # ---- Ashwell pathway ----
    ('rxn05673_c0', 'D-galacturonate proton symporter',
        {'2.A.1.14.-'}, ['galacturonate transporter', 'galacturonic acid transporter',
                         'gat1', 'gat', 'mfs', 'sugar transporter'], 'C. higginsianum'),
    ('rxn07491_c0', 'D-galacturonate reductase (NADPH; GAR1, fungal)',
        {'1.1.1.365'}, ['galacturonate reductase', 'gaaa', 'gar1', 'gar2'],
        'C. higginsianum'),
    ('rxn21749_c0', 'L-galactonate dehydratase (LGD1)',
        {'4.2.1.146'}, ['l-galactonate dehydratase', 'gaab', 'lgd1'],
        'C. higginsianum'),
    ('rxn21750_c0', '2-keto-3-deoxy-L-galactonate aldolase (LGA1)',
        {'4.1.2.55'}, ['2-keto-3-deoxy-l-galactonate aldolase', 'gaac', 'lga1',
                        '2-dehydro-3-deoxy-l-galactonate', 'kdg aldolase'],
        'C. higginsianum'),
    ('rxn09954_c0', 'L-glyceraldehyde reductase (GLD1)',
        {'1.1.1.21'}, ['glyceraldehyde reductase', 'gaad', 'gld1',
                        'aldose reductase', 'aldo-keto reductase'],
        'C. higginsianum'),
]

BASE      = '/home/janakae/fungalTemplate/imm904CobraModel'
GBFF      = f'{BASE}/gpr-update/C_higgensium.gbff'
BLAST_DB  = f'{BASE}/gpr-update/blast_db/fsp237'
OUT_DIR   = f'{BASE}/simulations/gapfill_v1_v2'
QUERY_FAA = f'{OUT_DIR}/candidates/v1_query_candidates.faa'
PROV_TSV  = f'{OUT_DIR}/candidates/v1_query_provenance.tsv'
BLAST_TSV = f'{OUT_DIR}/blast/v1_blast_results.tsv'
PICKS_TSV = f'{OUT_DIR}/blast/v1_rxn_to_gene_mapping.tsv'

PID_MIN  = 30.0
QCOV_MIN = 50.0
EVAL_MAX = 1e-10

os.makedirs(os.path.dirname(QUERY_FAA), exist_ok=True)
os.makedirs(os.path.dirname(BLAST_TSV), exist_ok=True)


def parse_gbff_proteins(path):
    """Yield {locus_tag, product, ec_numbers, protein_id, sequence} per CDS."""
    cur = None
    inside_seq = False
    in_cds = False
    cds_depth = 0
    in_translation = False
    with open(path) as fh:
        for line in fh:
            ls = line.rstrip()
            stripped = ls.strip()
            # New feature key indicator (column 6 has a token)
            if ls.startswith('     CDS '):
                if cur and cur.get('translation'):
                    yield cur
                cur = {'locus_tag': None, 'product': '', 'ec_numbers': set(),
                        'protein_id': None, 'translation': ''}
                in_cds = True
                in_translation = False
                continue
            if ls.startswith('     ') and not ls.startswith('                     '):
                # New non-CDS feature -> close out current CDS
                if cur and cur.get('translation'):
                    yield cur; cur = None
                in_cds = False
                in_translation = False
                continue
            if not cur:
                continue
            if not in_cds:
                continue
            # Qualifier lines start with 21 spaces and '/'
            if stripped.startswith('/'):
                # Close any in-progress multi-line value
                in_translation = False
                if stripped.startswith('/locus_tag='):
                    cur['locus_tag'] = stripped.split('=',1)[1].strip('"')
                elif stripped.startswith('/product='):
                    val = stripped.split('=',1)[1].strip('"')
                    cur['product'] = val
                    # Continued? leave as-is for now; product is usually 1-2 lines
                elif stripped.startswith('/EC_number='):
                    val = stripped.split('=',1)[1].strip('"')
                    cur['ec_numbers'].add(val)
                elif stripped.startswith('/protein_id='):
                    cur['protein_id'] = stripped.split('=',1)[1].strip('"')
                elif stripped.startswith('/translation='):
                    val = stripped.split('=',1)[1]
                    val = val.lstrip('"')
                    in_translation = True
                    # Could end on same line with closing quote
                    if val.endswith('"'):
                        cur['translation'] = val.rstrip('"')
                        in_translation = False
                    else:
                        cur['translation'] = val
            else:
                # Continuation
                if in_translation:
                    val = stripped
                    if val.endswith('"'):
                        cur['translation'] += val.rstrip('"')
                        in_translation = False
                    else:
                        cur['translation'] += val
                elif cur.get('product') and not cur['product'].endswith('"'):
                    # multi-line product -- just append a space and the text
                    # but only if the next line doesn't start a new qualifier
                    if not stripped.startswith('/'):
                        # Continuation if the value didn't close (still in product)
                        # Simplified -- we tolerate slightly truncated product strings
                        pass
        if cur and cur.get('translation'):
            yield cur


def score_protein(cds, ec_set, kw_list):
    """Higher = better match."""
    s = 0
    if cds['ec_numbers'] & ec_set:
        s += 100
    p = (cds['product'] or '').lower()
    for kw in kw_list:
        if kw.lower() in p:
            s += 10
    return s


def main():
    print(f'parsing GBFF: {GBFF}')
    cds_list = list(parse_gbff_proteins(GBFF))
    cds_list = [c for c in cds_list if c['translation']]
    print(f'  loaded {len(cds_list)} CDS records with protein seqs')

    # For each gapfill rxn, score and pick top-2 candidates
    queries = []  # list of (query_id, sequence, rxn_id, locus_tag, product, score, ec_set)
    print(f'\n=== selecting candidates ===')
    for rxn_id, label, ecs, kws, src in GAPFILL_TARGETS:
        scored = [(score_protein(c, ecs, kws), c) for c in cds_list]
        scored = [(s, c) for s, c in scored if s > 0]
        scored.sort(key=lambda x: -x[0])
        top = scored[:2]   # pick top 2 per reaction
        print(f'\n{rxn_id} ({label})')
        for s, c in top:
            qid = f'{rxn_id}__{c["locus_tag"]}'
            queries.append({
                'query_id': qid, 'sequence': c['translation'],
                'rxn_id': rxn_id, 'rxn_label': label,
                'locus_tag': c['locus_tag'], 'product': c['product'],
                'ec_numbers': ';'.join(sorted(c['ec_numbers'])),
                'protein_id': c['protein_id'] or '', 'source_strain': src,
                'score': s,
            })
            print(f'  + {c["locus_tag"]} (score {s}) "{c["product"][:55]}"  EC={sorted(c["ec_numbers"])}')
        if not top:
            print(f'  ! no C. higginsianum candidate matched')

    # Write query FASTA + provenance TSV
    print(f'\n=== writing {len(queries)} queries to {QUERY_FAA} ===')
    with open(QUERY_FAA, 'w') as fh:
        for q in queries:
            fh.write(f'>{q["query_id"]}\n')
            for line in textwrap.wrap(q['sequence'], 70):
                fh.write(line + '\n')
    import csv
    with open(PROV_TSV, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['query_id','rxn_id','rxn_label',
                                            'locus_tag','product','ec_numbers',
                                            'protein_id','source_strain','score'],
                            delimiter='\t')
        w.writeheader()
        for q in queries:
            row = {k: q.get(k, '') for k in w.fieldnames}
            w.writerow(row)

    # Run blastp
    print(f'\n=== running blastp ===')
    out = subprocess.run(
        ['blastp', '-query', QUERY_FAA, '-db', BLAST_DB,
         '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend '
                    'sstart send evalue bitscore qcovs qlen slen',
         '-max_target_seqs', '5', '-evalue', '1e-5'],
        capture_output=True, text=True, check=True)
    with open(BLAST_TSV, 'w') as fh:
        fh.write('qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend'
                 '\tsstart\tsend\tevalue\tbitscore\tqcovs\tqlen\tslen\n')
        fh.write(out.stdout)
    n_hits = len(out.stdout.strip().split('\n')) if out.stdout.strip() else 0
    print(f'  {n_hits} BLAST lines -> {BLAST_TSV}')

    # Parse + pick per rxn
    import collections
    queries_by_id = {q['query_id']: q for q in queries}
    hits_by_query = collections.defaultdict(list)
    for line in out.stdout.strip().split('\n'):
        if not line: continue
        parts = line.split('\t')
        d = dict(zip(['qseqid','sseqid','pident','length','mismatch','gapopen',
                       'qstart','qend','sstart','send','evalue','bitscore',
                       'qcovs','qlen','slen'], parts))
        d['pident'] = float(d['pident']); d['qcovs'] = float(d['qcovs'])
        d['evalue'] = float(d['evalue']); d['bitscore'] = float(d['bitscore'])
        hits_by_query[d['qseqid']].append(d)

    # Build per-rxn pick: best gene_* hit that passes thresholds across all queries for that rxn
    picks = {}  # rxn_id -> dict
    for q in queries:
        rxn = q['rxn_id']
        for h in hits_by_query.get(q['query_id'], []):
            if not h['sseqid'].startswith('gene_'):  # skip MSTRG
                continue
            passes = (h['pident'] >= PID_MIN and h['qcovs'] >= QCOV_MIN
                      and h['evalue'] <= EVAL_MAX)
            if not passes: continue
            cur = picks.get(rxn)
            if cur is None or h['bitscore'] > cur['bitscore']:
                picks[rxn] = {
                    'rxn_id': rxn, 'rxn_label': q['rxn_label'],
                    'source_strain': q['source_strain'],
                    'source_locus_tag': q['locus_tag'],
                    'source_product': q['product'],
                    'source_ec': q['ec_numbers'],
                    'source_protein_id': q['protein_id'],
                    'fsp237_gene': h['sseqid'],
                    'pident': h['pident'], 'qcovs': h['qcovs'],
                    'evalue': h['evalue'], 'bitscore': h['bitscore'],
                    'length': h['length'], 'qlen': h['qlen'], 'slen': h['slen'],
                }

    # Confidence classification
    for rxn in picks:
        p = picks[rxn]
        if p['pident'] >= 60 and p['qcovs'] >= 70:
            p['confidence'] = 'high'
        elif p['pident'] >= 45 and p['qcovs'] >= 60:
            p['confidence'] = 'medium'
        else:
            p['confidence'] = 'low'
        p['note'] = ''

    # Reactions with no pick
    for rxn_id, label, _, _, src in GAPFILL_TARGETS:
        if rxn_id not in picks:
            picks[rxn_id] = {
                'rxn_id': rxn_id, 'rxn_label': label,
                'source_strain': src, 'source_locus_tag': '',
                'source_product': '', 'source_ec': '',
                'source_protein_id': '',
                'fsp237_gene': '', 'pident': '', 'qcovs': '',
                'evalue': '', 'bitscore': '', 'length': '', 'qlen': '', 'slen': '',
                'confidence': 'none', 'note': 'No C. higg candidate or no BLAST hit passes thresholds',
            }

    # Write final mapping TSV
    with open(PICKS_TSV, 'w', newline='') as fh:
        cols = ['rxn_id','rxn_label','source_strain','source_locus_tag','source_product',
                'source_ec','source_protein_id','fsp237_gene','pident','qcovs',
                'evalue','bitscore','length','qlen','slen','confidence','note']
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        w.writeheader()
        for rxn_id, _, _, _, _ in GAPFILL_TARGETS:
            w.writerow({k: picks[rxn_id].get(k, '') for k in cols})
    print(f'\n=== final mapping ===')
    for rxn_id, label, _, _, _ in GAPFILL_TARGETS:
        p = picks[rxn_id]
        gene = p['fsp237_gene'] or 'NONE'
        pid  = f'{p["pident"]}' if p['pident'] != '' else '-'
        qcov = f'{p["qcovs"]}'   if p['qcovs']  != '' else '-'
        print(f'  {rxn_id:30s} -> {gene:15s}  pid={pid:>5}  qcov={qcov:>5}  conf={p["confidence"]}')
    print(f'\nsaved: {PICKS_TSV}')


if __name__ == '__main__':
    main()
