#!/opt/env/modelseed/bin/python3
"""Assemble the compact JSON payload the site's RNA-seq tab consumes.

Emits ROOT/outputs/rnaseq_payload.json with keys:
  coverage        — Stage-0 coverage stats
  condition_fit   — Stage-1 per-condition summary (aerobic + anaerobic pivoted)
  per_reaction    — Stage-0 reaction expression scores + Stage-1 default-condition flags
  stage2_summary  — GIMME/iMAT/E-Flux per (cond, O2, method, threshold)
  orphan_priority — Stage-4 cross-condition orphan reactions ranked
  thresholds      — the exact quantile cutoffs used
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path('/home/janakae/fungalTemplate/imm904CobraModel/rna_seq_integration')
OUTPUTS = ROOT / 'outputs'
PAYLOAD_JSON = OUTPUTS / 'rnaseq_payload.json'

DEFAULT_CONDITION_FOR_PER_REACTION = '06_glucose_low'  # user's biological anchor


def _pivot_conditions(df, method_filter=None, threshold_filter=None):
    """Pivot per-(cond,O2) rows into per-condition rows with aer+ana columns."""
    if method_filter is not None:
        df = df[df['method'] == method_filter] if 'method' in df.columns else df
    if threshold_filter is not None and 'threshold' in df.columns:
        df = df[df['threshold'] == threshold_filter]
    by_cond = {}
    order = []
    for _, row in df.iterrows():
        cid = row['condition_id']
        if cid not in by_cond:
            by_cond[cid] = {
                'condition_id': cid,
                'label': row.get('label', cid),
                'stage': row.get('stage', ''),
            }
            order.append(cid)
        o2 = row['O2']
        for k in row.index:
            if k in ('condition_id', 'label', 'stage', 'O2', 'notes',
                     'method', 'threshold'):
                continue
            val = row[k]
            if pd.isna(val):
                val = None
            elif hasattr(val, 'item'):
                val = val.item()
            by_cond[cid][f'{k}_{o2[:3]}'] = val
    return [by_cond[c] for c in order]


def build():
    coverage = json.loads((OUTPUTS / 'coverage_summary.json').read_text())

    stage1_summary = pd.read_csv(OUTPUTS / 'stage1_summary.tsv', sep='\t')
    condition_fit = _pivot_conditions(stage1_summary)

    # Per-reaction table: expression scores + default-condition Stage 1 flags
    rexp = pd.read_csv(OUTPUTS / 'reaction_expression.tsv', sep='\t', comment='#')
    st1_default_path = OUTPUTS / 'stage1_overlay' / f'{DEFAULT_CONDITION_FOR_PER_REACTION}_aerobic.tsv'
    if st1_default_path.exists():
        st1_default = pd.read_csv(st1_default_path, sep='\t')[
            ['rxn_id', 'flux_pFBA', 'agreement']
        ].rename(columns={'flux_pFBA': 'flux_glucose_low_aer',
                           'agreement': 'agreement_glucose_low_aer'})
        rexp = rexp.merge(st1_default, on='rxn_id', how='left')
    per_reaction = []
    for _, r in rexp.iterrows():
        row = {
            'rxn_id': r['rxn_id'],
            'name': r.get('name', '') or '',
            'subsystem': r.get('subsystem', '') or '',
            'gpr': r.get('gpr', '') or '',
            'n_genes': int(r['n_genes']) if pd.notna(r['n_genes']) else 0,
            'n_genes_with_expr': int(r['n_genes_with_expr']) if pd.notna(r['n_genes_with_expr']) else 0,
            'agg_mean_log2TPMp1': None if pd.isna(r.get('agg_mean_log2TPMp1')) else round(float(r['agg_mean_log2TPMp1']), 3),
            'expression_bin': r.get('expression_bin', 'absent'),
            'notes': r.get('notes', '') or '',
        }
        if 'flux_glucose_low_aer' in r:
            row['flux_glucose_low_aer'] = (None if pd.isna(r['flux_glucose_low_aer'])
                                            else round(float(r['flux_glucose_low_aer']), 4))
            row['agreement_glucose_low_aer'] = r.get('agreement_glucose_low_aer', '')
        per_reaction.append(row)

    # Stage 2 summary — pivot method+threshold into rows
    stage2 = None
    stage2_path = OUTPUTS / 'stage2_summary.tsv'
    if stage2_path.exists():
        stage2_df = pd.read_csv(stage2_path, sep='\t')
        stage2 = []
        for _, r in stage2_df.iterrows():
            stage2.append({
                'condition_id': r['condition_id'],
                'label': r.get('label', r['condition_id']),
                'stage': r.get('stage', ''),
                'O2': r['O2'], 'method': r['method'],
                'threshold': r.get('threshold', ''),
                'status': r['status'],
                'biomass': None if pd.isna(r['biomass']) else round(float(r['biomass']), 4),
                'n_active': int(r['n_active']) if pd.notna(r['n_active']) else 0,
                'n_H_active': None if pd.isna(r.get('n_H_active')) else int(r['n_H_active']),
                'n_H_total': None if pd.isna(r.get('n_H_total')) else int(r['n_H_total']),
                'n_L_inactive': None if pd.isna(r.get('n_L_inactive')) else int(r['n_L_inactive']),
                'n_L_total': None if pd.isna(r.get('n_L_total')) else int(r['n_L_total']),
            })

    # Stage 4 orphan priority
    orphan_priority = []
    op_path = OUTPUTS / 'orphan_priority.tsv'
    if op_path.exists():
        opdf = pd.read_csv(op_path, sep='\t')
        for _, r in opdf.head(50).iterrows():
            orphan_priority.append({
                'rxn_id': r['rxn_id'],
                'name': r.get('name', '') or '',
                'subsystem': r.get('subsystem', '') or '',
                'equation': r.get('equation', '') or '',
                'n_conditions_active': int(r['n_conditions_active']),
                'max_abs_flux': round(float(r['max_abs_flux']), 4),
                'sum_abs_flux': round(float(r['sum_abs_flux']), 4),
                'top_candidates': r.get('top_candidates', '') or '',
                'top_candidate_expr': r.get('top_candidate_expr', '') or '',
            })

    payload = {
        'source': 'rna_seq_integration/',
        'dataset_label': 'S1 (3 replicates, G1/G2/G3)',
        'note': 'Single-condition RNA-seq (S1). Cross-condition MADE deferred until additional biological groups arrive.',
        'coverage': coverage,
        'default_condition_for_per_reaction': DEFAULT_CONDITION_FOR_PER_REACTION,
        'condition_fit': condition_fit,
        'per_reaction': per_reaction,
        'stage2_summary': stage2,
        'orphan_priority': orphan_priority,
    }
    PAYLOAD_JSON.write_text(json.dumps(payload))
    print(f'wrote {PAYLOAD_JSON}')
    print(f'  coverage           : {coverage["model_genes_with_expression"]} / {coverage["model_genes"]} model genes')
    print(f'  condition_fit rows : {len(condition_fit)}')
    print(f'  per_reaction rows  : {len(per_reaction)}')
    print(f'  stage2 rows        : {len(stage2) if stage2 else 0}')
    print(f'  orphan rows        : {len(orphan_priority)} (top-50)')


if __name__ == '__main__':
    build()
