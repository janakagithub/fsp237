#!/usr/bin/env python3
"""COMPREHENSIVE GPR update for the FSP237 model.

Combines:
  - CH63R BLAST mapping (from both Excels)
  - Excel curated gene sets per reaction (Integration_Higgginsianum_GEM - 4JK.xlsx
    + 15June2026_Missing reactions GEM_higginsianum.xlsx)

For EVERY reaction in the model:
  1. Keep all `gene_*` tokens (FSP237 native; user instruction).
  2. Look up the reaction in the Excel sheets by `rxn id` and get its
     curated CH63R gene set.
  3. Map each CH63R via the BLAST table to a gene_* (above quality threshold);
     dedupe.
  4. Remove all yeast-style (Y####W/C) and accession-style (NP_*) tokens that
     have no curated CH63R equivalent for this reaction.
  5. Build new GPR = gene_* (preserved) + mapped(CH63R via Excel) -- "or" join.
  6. Save the model + change log + summary.
"""
import json
import os
import re
import sys

import cobra
import pandas as pd

GPR_DIR = '/home/janakae/fungalTemplate/imm904CobraModel/gpr-update'
INPUT_MODEL = '/home/janakae/fungalTemplate/imm904CobraModel/fsp237_biomass_extension/fsp237_atp_safe_gsm_extended.json'
OUTPUT_MODEL = f'{GPR_DIR}/fsp237_atp_safe_gsm_gpr_updated.json'
BLAST_TSV = f'{GPR_DIR}/blast_results_all.tsv'
MAPPING_TSV = f'{GPR_DIR}/ch63r_to_fsp237_mapping_full.tsv'
CHANGELOG = f'{GPR_DIR}/gpr_change_log.tsv'
FLAGGED = f'{GPR_DIR}/flagged_reactions.tsv'

XLS = [
    '/home/janakae/fungalTemplate/imm904CobraModel/Integration_Higgginsianum_GEM - 4JK.xlsx',
    '/home/janakae/fungalTemplate/imm904CobraModel/fsp237_biomass_extension/15June2026_Missing reactions GEM_higginsianum.xlsx',
]

# Quality thresholds for "this CH63R has a reliable FSP237 ortholog"
PID_THRESH = 30.0
QCOV_THRESH = 50.0
EVAL_THRESH = 1e-10

# Token regex
GENE_RE = re.compile(r'^gene_\d+$')
CH_RE = re.compile(r'^CH63R_\d+$')


# ---- 1. Build the BLAST mapping (CH63R -> best gene_*) ---------------------

def build_blast_mapping():
    cols = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend',
            'sstart','send','evalue','bitscore','qcovs','qlen','slen']
    df = pd.read_csv(BLAST_TSV, sep='\t', names=cols, header=None)
    df['is_gene'] = df['sseqid'].str.startswith('gene_')
    df['is_mstrg'] = df['sseqid'].str.contains('MSTRG')
    df['passes'] = ((df['pident'] >= PID_THRESH) & (df['qcovs'] >= QCOV_THRESH)
                    & (df['evalue'] <= EVAL_THRESH))

    rows = []
    mapping = {}
    weak = set()
    for q in df['qseqid'].unique():
        sub = df[df['qseqid'] == q].sort_values('bitscore', ascending=False)
        best = sub.iloc[0]
        gene_pass = sub[sub['is_gene'] & sub['passes']].head(1)
        mstrg_best = sub[sub['is_mstrg']].head(1)
        mapped = gene_pass['sseqid'].iloc[0] if len(gene_pass) else None
        flag = ''
        if mapped is None:
            flag = 'NO_GENE_HIT_PASSES_THRESHOLDS'
        elif gene_pass['pident'].iloc[0] < 50:
            flag = 'WEAK_MATCH (pident < 50%)'
            weak.add(q)
        if mapped:
            mapping[q] = mapped
        rows.append({
            'ch63r_id': q,
            'mapped_gene': mapped,
            'mapped_pident': gene_pass['pident'].iloc[0] if len(gene_pass) else None,
            'mapped_qcovs': gene_pass['qcovs'].iloc[0] if len(gene_pass) else None,
            'mapped_evalue': gene_pass['evalue'].iloc[0] if len(gene_pass) else None,
            'mapped_bitscore': gene_pass['bitscore'].iloc[0] if len(gene_pass) else None,
            'best_overall_hit': best['sseqid'],
            'best_overall_pident': best['pident'],
            'best_overall_is_mstrg': bool(best['is_mstrg']),
            'mstrg_best_hit': mstrg_best['sseqid'].iloc[0] if len(mstrg_best) else None,
            'mstrg_best_pident': mstrg_best['pident'].iloc[0] if len(mstrg_best) else None,
            'flag': flag,
        })
    pd.DataFrame(rows).to_csv(MAPPING_TSV, sep='\t', index=False)
    print(f'BLAST mapping: {len(mapping)} CH63R->gene_* / {len(rows)} CH63R queries')
    print(f'  flagged weak    : {len(weak)}')
    print(f'  no gene hit     : {sum(1 for r in rows if r["flag"]=="NO_GENE_HIT_PASSES_THRESHOLDS")}')
    return mapping, weak


# ---- 2. Build per-reaction curated CH63R set from Excel sheets -------------

def build_rxn_ch63r_index():
    """For each rxn id in the Excels, what CH63R genes were curated?"""
    out = {}
    for xl in XLS:
        xf = pd.ExcelFile(xl)
        for sh in xf.sheet_names:
            df = pd.read_excel(xl, sheet_name=sh, dtype=str).fillna('')
            if 'rxn id' not in df.columns: continue
            gene_col = next((c for c in df.columns if 'gene ID' in c or 'gene id' in c.lower()), None)
            if not gene_col: continue
            for _, row in df.iterrows():
                rid = str(row['rxn id']).strip()
                if not rid or rid == 'nan' or rid == 'rxn': continue
                gene_blob = str(row[gene_col])
                chs = re.findall(r'CH63R_\d+', gene_blob)
                if not chs: continue
                # Append to set (may merge across sheets)
                out.setdefault(rid, set()).update(chs)
    print(f'Excel curated rxns with CH63R genes: {len(out)}')
    return out


# ---- 3. Update each reaction's GPR -----------------------------------------

def update_model_gprs(model, blast_map, weak, rxn_ch63r):
    changelog = []
    flagged = []

    for r in model.reactions:
        gpr = r.gene_reaction_rule or ''
        # 1. Extract gene_* tokens to preserve (FSP237 native -- never lose these)
        tokens = re.findall(r'[A-Za-z][A-Za-z0-9_.]*', gpr)
        tokens = [t for t in tokens if t not in ('and', 'or')]
        kept_gene = [t for t in tokens if GENE_RE.match(t)]
        foreign = [t for t in tokens if not GENE_RE.match(t)]

        # 2. Pull curated CH63R for this reaction (if Excel has an entry)
        curated_ch = rxn_ch63r.get(r.id, set())
        # also pick up CH63R tokens that may have been in the GPR itself
        for t in foreign:
            if CH_RE.match(t): curated_ch.add(t)
        # map curated CH63R -> gene_*
        mapped = []
        unmapped = []
        weak_used = []
        for ch in sorted(curated_ch):
            if ch in blast_map:
                mapped.append(blast_map[ch])
                if ch in weak: weak_used.append(ch)
            else:
                unmapped.append(ch)

        # 3. Dedupe + assemble new GPR
        new_gene_set = list(dict.fromkeys(kept_gene + mapped))
        new_gpr = ' or '.join(new_gene_set) if new_gene_set else ''

        # 4. Only update if it changed
        if new_gpr == gpr: continue

        # 5. Record change
        had_foreign = bool(foreign)
        n_removed_yeast = sum(1 for t in foreign if not CH_RE.match(t))
        changelog.append({
            'rxn_id': r.id, 'name': r.name,
            'old_gpr': gpr, 'new_gpr': new_gpr,
            'kept_gene_count': len(kept_gene),
            'mapped_from_ch63r': len(mapped),
            'unmapped_ch63r_removed': len(unmapped),
            'foreign_tokens_removed': n_removed_yeast,
            'curated_ch63r_from_excel': ';'.join(sorted(curated_ch)),
            'foreign_tokens_in_gpr': ';'.join(t for t in foreign if not CH_RE.match(t)),
            'weak_matches_used': ';'.join(weak_used),
        })

        # 6. Flag if needed
        msg = []
        if unmapped:
            msg.append(f'unmapped CH63R dropped: {sorted(unmapped)}')
        if weak_used:
            msg.append(f'weak BLAST matches used: {weak_used}')
        if not new_gpr and gpr:
            msg.append('GPR became EMPTY')
        if msg:
            flagged.append({
                'rxn_id': r.id, 'name': r.name,
                'old_gpr': gpr, 'new_gpr': new_gpr,
                'issue': ' | '.join(msg),
            })

        r.gene_reaction_rule = new_gpr

    return changelog, flagged


def main():
    print(f'loading: {INPUT_MODEL}')
    m = cobra.io.load_json_model(INPUT_MODEL)
    print(f'  start: {len(m.reactions)} reactions, {len(m.genes)} genes')

    print('\n=== 1. building BLAST mapping ===')
    blast_map, weak = build_blast_mapping()

    print('\n=== 2. building per-reaction CH63R index from Excel ===')
    rxn_ch63r = build_rxn_ch63r_index()

    print('\n=== 3. updating model GPRs ===')
    changelog, flagged = update_model_gprs(m, blast_map, weak, rxn_ch63r)

    # Sanity: count remaining foreign tokens
    n_foreign_after = 0
    n_rxn_with_foreign_after = 0
    for r in m.reactions:
        if not r.gene_reaction_rule: continue
        toks = [t for t in re.findall(r'[A-Za-z][A-Za-z0-9_.]*', r.gene_reaction_rule)
                if t not in ('and','or')]
        bad = [t for t in toks if not GENE_RE.match(t)]
        if bad:
            n_foreign_after += len(bad)
            n_rxn_with_foreign_after += 1

    cobra.io.save_json_model(m, OUTPUT_MODEL)
    pd.DataFrame(changelog).to_csv(CHANGELOG, sep='\t', index=False)
    pd.DataFrame(flagged).to_csv(FLAGGED, sep='\t', index=False)

    print('\n=== SUMMARY ===')
    print(f'  reactions updated     : {len(changelog)}')
    print(f'  reactions flagged     : {len(flagged)}')
    print(f'  total CH63R mapped    : {sum(c["mapped_from_ch63r"] for c in changelog)}')
    print(f'  total foreign removed : {sum(c["foreign_tokens_removed"] for c in changelog)}')
    print(f'  total unmapped dropped: {sum(c["unmapped_ch63r_removed"] for c in changelog)}')
    print(f'\n  REMAINING foreign tokens in any GPR: {n_foreign_after}')
    print(f'  REMAINING reactions with foreign tokens: {n_rxn_with_foreign_after}')

    # Smoke test
    sol = m.optimize()
    print(f'\n  model optimizes to: {sol.objective_value:.6g}')

    print(f'\nsaved: {OUTPUT_MODEL}')
    print(f'saved: {CHANGELOG}')
    print(f'saved: {FLAGGED}')

if __name__ == '__main__':
    main()
