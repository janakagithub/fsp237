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
import math
from pathlib import Path

import pandas as pd


def _clean_nan(obj):
    """Recursively replace NaN / +/-Inf floats with None so json.dump
    produces browser-parseable JSON (JavaScript's JSON.parse rejects bare
    NaN and Infinity)."""
    if isinstance(obj, dict):
        return {k: _clean_nan(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_clean_nan(v) for v in obj]
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj

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

    # Stage 5 pathway analysis
    pathway_summary = []
    ps_path = OUTPUTS / 'pathway_summary.tsv'
    if ps_path.exists():
        for _, r in pd.read_csv(ps_path, sep='\t').iterrows():
            pathway_summary.append({
                'pathway': r['pathway'],
                'n_reactions': int(r['n_reactions']),
                'n_with_gpr': int(r['n_with_gpr']),
                'n_scored': int(r['n_scored']),
                'mean_expression': None if pd.isna(r['mean_expression']) else round(float(r['mean_expression']), 3),
                'median_expression': None if pd.isna(r['median_expression']) else round(float(r['median_expression']), 3),
                'hi_n': int(r['hi_n']),
                'med_n': int(r['med_n']),
                'lo_n': int(r['lo_n']),
                'absent_n': int(r['absent_n']),
                'hi_frac': round(float(r['hi_frac']), 3),
                'odds_ratio': round(float(r['odds_ratio']), 3),
                'p_fisher_greater': round(float(r['p_fisher_greater']), 5),
                'bh_q': round(float(r['bh_q']), 4),
            })

    pathway_matrix = []
    pm_path = OUTPUTS / 'pathway_condition_matrix.tsv'
    if pm_path.exists():
        for _, r in pd.read_csv(pm_path, sep='\t').iterrows():
            pathway_matrix.append({
                'pathway': r['pathway'],
                'condition_id': r['condition_id'],
                'stage': r['stage'],
                'O2': r['O2'],
                'sum_abs_flux': round(float(r['sum_abs_flux']), 3),
                'n_active': int(r['n_active']),
                'n_total': int(r['n_total']),
                'active_rate': round(float(r['active_rate']), 3),
            })

    reporter_metabolites = []
    rm_path = OUTPUTS / 'reporter_metabolites.tsv'
    if rm_path.exists():
        for _, r in pd.read_csv(rm_path, sep='\t').head(60).iterrows():
            reporter_metabolites.append({
                'metabolite_id': r['metabolite_id'],
                'name': r.get('name', '') or '',
                'compartment': r.get('compartment', '') or '',
                'formula': r.get('formula', '') or '',
                'k_reactions': int(r['k_reactions']),
                'z_reporter': round(float(r['z_reporter']), 3),
                'p_one_sided': round(float(r['p_one_sided']), 4),
            })

    biomass_corr = []
    bc_path = OUTPUTS / 'biomass_pathway_corr.tsv'
    if bc_path.exists():
        for _, r in pd.read_csv(bc_path, sep='\t').iterrows():
            biomass_corr.append({
                'pathway': r['pathway'],
                'n_conditions': int(r['n_conditions']),
                'pearson_r_flux_vs_biomass': round(float(r['pearson_r_flux_vs_biomass']), 3),
                'p_value': round(float(r['p_value']), 4),
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
        'pathway_summary': pathway_summary,
        'pathway_matrix': pathway_matrix,
        'reporter_metabolites': reporter_metabolites,
        'biomass_pathway_corr': biomass_corr,
    }
    payload = _clean_nan(payload)
    # allow_nan=False makes json.dump raise if any NaN/Inf slipped through
    # the sanitiser instead of silently emitting invalid JSON.
    PAYLOAD_JSON.write_text(json.dumps(payload, allow_nan=False))
    print(f'wrote {PAYLOAD_JSON}')
    print(f'  coverage           : {coverage["model_genes_with_expression"]} / {coverage["model_genes"]} model genes')
    print(f'  condition_fit rows : {len(condition_fit)}')
    print(f'  per_reaction rows  : {len(per_reaction)}')
    print(f'  stage2 rows        : {len(stage2) if stage2 else 0}')
    print(f'  orphan rows        : {len(orphan_priority)} (top-50)')


if __name__ == '__main__':
    build()
