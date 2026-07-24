#!/opt/env/modelseed/bin/python3
"""Stage 1 — pFBA × expression overlay per condition × O2.

For each (condition, O2) in the 18-condition panel:
  - reapply the exact same medium as run_simulation_panel.py
  - solve pFBA (parsimonious FBA, cobra.flux_analysis.pfba)
  - merge with reaction_expression.tsv (Stage 0)
  - classify each reaction into an agreement category:
        SUPPORTED                    flux != 0 & expr_bin in {hi, med}
        WEAK_SUPPORT                 flux != 0 & expr_bin == 'lo'
        CONFLICT_FLUX_NO_EXPR        flux != 0 & expr_bin == 'absent' & has GPR
        ORPHAN_FLUX                  flux != 0 & no GPR
        PRIMED_NOT_USED              flux == 0 & expr_bin == 'hi'
        SILENT_OK                    everything else with flux == 0

Emit:
  - outputs/stage1_overlay/<cond>_<o2>.tsv  (per-reaction)
  - outputs/stage1_summary.tsv              (per-condition roll-up)

The condition-fit `agreement_score` (higher = the S1 transcriptome fits that
medium better) is:
     SUPPORTED / (SUPPORTED + CONFLICT_FLUX_NO_EXPR + PRIMED_NOT_USED)
"""
import sys
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from cobra.flux_analysis import pfba

sys.path.insert(0, '/home/janakae/fsp237/simulations')
from run_simulation_panel import (  # noqa: E402
    CONDITIONS, INORGANIC, O2, AA_EXCHANGES, apply_media,
)

ROOT = Path('/home/janakae/fungalTemplate/imm904CobraModel/rna_seq_integration')
MODEL_PATH = Path('/home/janakae/fsp237/simulations/gapfill_v1_v2/models/'
                  'fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json')
REXP_TSV = ROOT / 'outputs' / 'reaction_expression.tsv'
OUT_DIR = ROOT / 'outputs' / 'stage1_overlay'
SUMMARY = ROOT / 'outputs' / 'stage1_summary.tsv'

FLUX_EPS = 1e-6


def classify(flux, bin_, has_gpr):
    nonzero = abs(flux) > FLUX_EPS
    if nonzero:
        if not has_gpr:
            return 'ORPHAN_FLUX'
        if bin_ in ('hi', 'med'):
            return 'SUPPORTED'
        if bin_ == 'lo':
            return 'WEAK_SUPPORT'
        return 'CONFLICT_FLUX_NO_EXPR'
    if bin_ == 'hi':
        return 'PRIMED_NOT_USED'
    return 'SILENT_OK'


def main():
    print(f'loading model : {MODEL_PATH.name}')
    model = cobra.io.load_json_model(str(MODEL_PATH))

    rexp = pd.read_csv(REXP_TSV, sep='\t', comment='#').set_index('rxn_id')
    print(f'loaded reaction_expression: {len(rexp)} rows')

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for cond_id, label, stage, c_sources, notes in CONDITIONS:
        for aerobic in (True, False):
            o2_tag = 'aerobic' if aerobic else 'anaerobic'
            with model:
                apply_media(model, c_sources, aerobic)
                for r in model.reactions:
                    r.objective_coefficient = 1 if r.id == 'bio_gsm' else 0
                sol = model.optimize()
                status = sol.status
                biomass = sol.objective_value if status == 'optimal' else None

                if biomass is None or biomass < FLUX_EPS:
                    print(f'  {cond_id:24s} {o2_tag[:3].upper():>3s}: '
                          f'infeasible / no growth   status={status}   bio={biomass}')
                    summary_rows.append({
                        'condition_id': cond_id, 'label': label, 'stage': stage,
                        'O2': o2_tag, 'status': status, 'biomass': biomass or 0.0,
                        'n_flux_nonzero': 0, 'n_gpr': int((rexp['n_genes'] > 0).sum()),
                        'SUPPORTED': 0, 'WEAK_SUPPORT': 0,
                        'CONFLICT_FLUX_NO_EXPR': 0, 'ORPHAN_FLUX': 0,
                        'PRIMED_NOT_USED': 0, 'SILENT_OK': 0,
                        'agreement_score': 0.0,
                    })
                    continue

                try:
                    psol = pfba(model)
                    fluxes = psol.fluxes
                except Exception as e:
                    print(f'  {cond_id:24s} {o2_tag[:3].upper():>3s}: '
                          f'pFBA failed ({e}), falling back to plain FBA')
                    fluxes = sol.fluxes

            # Build per-reaction table
            per = rexp.copy()
            per['flux_pFBA'] = per.index.map(lambda rid: float(fluxes.get(rid, 0.0)))
            per['has_gpr'] = per['n_genes'] > 0
            per['agreement'] = [
                classify(f, b, hg)
                for f, b, hg in zip(per['flux_pFBA'], per['expression_bin'], per['has_gpr'])
            ]
            per_reset = per.reset_index()
            per_out = OUT_DIR / f'{cond_id}_{o2_tag}.tsv'
            per_reset.to_csv(per_out, sep='\t', index=False)

            counts = per['agreement'].value_counts().to_dict()
            supported = counts.get('SUPPORTED', 0)
            conflict = counts.get('CONFLICT_FLUX_NO_EXPR', 0)
            primed = counts.get('PRIMED_NOT_USED', 0)
            denom = supported + conflict + primed
            score = supported / denom if denom else 0.0
            n_nonzero = int((per['flux_pFBA'].abs() > FLUX_EPS).sum())

            # ---- additional condition-fit metrics ----
            gpr_mask = per['n_genes'] > 0
            expr_scores = per.loc[gpr_mask, 'agg_mean_log2TPMp1'].to_numpy(dtype=float)
            abs_flux = per.loc[gpr_mask, 'flux_pFBA'].abs().to_numpy(dtype=float)
            expr_valid = ~np.isnan(expr_scores)
            if expr_valid.sum() > 10 and abs_flux[expr_valid].max() > FLUX_EPS:
                rho, _ = spearmanr(expr_scores[expr_valid], abs_flux[expr_valid])
                if np.isnan(rho):
                    rho = 0.0
            else:
                rho = 0.0

            hi_rxns = per[per['expression_bin'] == 'hi']
            n_hi = len(hi_rxns)
            n_hi_active = int((hi_rxns['flux_pFBA'].abs() > FLUX_EPS).sum())
            hi_recall = n_hi_active / n_hi if n_hi else 0.0

            gpr_active = per[gpr_mask & (per['flux_pFBA'].abs() > FLUX_EPS)]
            n_active_gpr = len(gpr_active)
            n_active_hi = int((gpr_active['expression_bin'].isin(['hi', 'med'])).sum())
            flux_precision = n_active_hi / n_active_gpr if n_active_gpr else 0.0

            summary_rows.append({
                'condition_id': cond_id, 'label': label, 'stage': stage,
                'O2': o2_tag, 'status': status, 'biomass': round(biomass, 6),
                'n_flux_nonzero': n_nonzero,
                'n_gpr': int((per['n_genes'] > 0).sum()),
                'SUPPORTED': supported,
                'WEAK_SUPPORT': counts.get('WEAK_SUPPORT', 0),
                'CONFLICT_FLUX_NO_EXPR': conflict,
                'ORPHAN_FLUX': counts.get('ORPHAN_FLUX', 0),
                'PRIMED_NOT_USED': primed,
                'SILENT_OK': counts.get('SILENT_OK', 0),
                'agreement_score': round(score, 4),
                'spearman_expr_vs_flux': round(float(rho), 4),
                'hi_expr_recall': round(hi_recall, 4),
                'flux_precision_hi_med': round(flux_precision, 4),
            })
            print(f'  {cond_id:24s} {o2_tag[:3].upper():>3s}: '
                  f'bio={biomass:.4f}  nonzero={n_nonzero:>4d}  '
                  f'SUP={supported:>4d} PNU={primed:>3d}  '
                  f'ρ={rho:+.3f}  hi_recall={hi_recall:.3f}  prec={flux_precision:.3f}')

    sdf = pd.DataFrame(summary_rows)
    sdf = sdf.sort_values(['O2', 'agreement_score'], ascending=[True, False])
    sdf.to_csv(SUMMARY, sep='\t', index=False)
    print(f'\nwrote {SUMMARY}')

    print('\n=== AEROBIC RANKING (biological expectation: 06_glucose_low high) ===')
    aer = sdf[sdf.O2 == 'aerobic'].copy()
    for col in ('agreement_score', 'spearman_expr_vs_flux', 'hi_expr_recall',
                'flux_precision_hi_med'):
        top = aer.sort_values(col, ascending=False).head(5)
        print(f'\n  by {col}:')
        for _, row in top.iterrows():
            marker = '  <-- glucose_low' if row.condition_id == '06_glucose_low' else ''
            print(f'    {row.condition_id:26s} {row.stage:16s} {row[col]:+.4f}{marker}')
    print('\n  06_glucose_low position across metrics:')
    for col in ('agreement_score', 'spearman_expr_vs_flux', 'hi_expr_recall',
                'flux_precision_hi_med'):
        aer_s = aer.sort_values(col, ascending=False).reset_index(drop=True)
        rank = int(aer_s.index[aer_s.condition_id == '06_glucose_low'][0]) + 1
        print(f'    {col:26s} rank {rank:>2d} / {len(aer_s)}')


if __name__ == '__main__':
    main()
