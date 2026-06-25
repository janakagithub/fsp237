#!/usr/bin/env python3
"""Re-run the 18-condition simulation panel against a chosen model and
record the results.

Usage:
  python test_gapfilled_model.py <model_json> <out_suffix>

e.g.
  python test_gapfilled_model.py models/fsp237_gapfilled_Version1_noGenes.json v1
  python test_gapfilled_model.py models/fsp237_gapfilled_Version2_gapfill_genes_integrated.json v2

Outputs (relative to this directory):
  reports/<suffix>_simulation_results.tsv
  reports/<suffix>_per_condition/<cond_id>_aerobic.tsv  etc.
"""
import os
import sys

import cobra
import pandas as pd

# Reuse the conditions + helpers from the panel driver
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))   # add simulations/ to path
from run_simulation_panel import (
    CONDITIONS, INORGANIC, O2, AA_EXCHANGES,
    apply_media, carbon_atoms_per_cpd, total_C_uptake,
    actual_C_uptake_from_solution,
)


def run_panel(model_path, suffix):
    out_dir  = f'{HERE}/reports'
    per_dir  = f'{out_dir}/{suffix}_per_condition'
    summary  = f'{out_dir}/{suffix}_simulation_results.tsv'
    os.makedirs(per_dir, exist_ok=True)

    print(f'loading: {model_path}')
    model = cobra.io.load_json_model(model_path)
    print(f'  {len(model.reactions)} rxns / {len(model.metabolites)} mets / {len(model.genes)} genes')
    c_atoms = carbon_atoms_per_cpd(model)

    rows = []
    print(f'\nrunning {len(CONDITIONS)} conditions x 2 (aerobic/anaerobic)')
    for cond_id, label, stage, c_sources, notes in CONDITIONS:
        for aerobic in (True, False):
            with model:
                missing = apply_media(model, c_sources, aerobic)
                for r in model.reactions:
                    r.objective_coefficient = 1 if r.id == 'bio_gsm' else 0
                sol = model.optimize()
                biomass = sol.objective_value if sol.status == 'optimal' else None
                status  = sol.status
                active = sum(1 for v in sol.fluxes.values if abs(v) > 1e-6) if biomass else 0
                c_in_theor = total_C_uptake(model, c_sources, c_atoms)
                c_in_real  = actual_C_uptake_from_solution(sol, c_sources, c_atoms) if biomass else 0.0
                bio_per_C  = (biomass / c_in_real) if (biomass and c_in_real > 0) else 0.0
                if biomass and biomass > 1e-6:
                    df = pd.DataFrame({'rxn_id': sol.fluxes.index, 'flux': sol.fluxes.values})
                    df = df[df['flux'].abs() > 1e-6].sort_values('flux', key=lambda c: c.abs(),
                                                                  ascending=False)
                    df.to_csv(f'{per_dir}/{cond_id}_{"aerobic" if aerobic else "anaerobic"}.tsv',
                              sep='\t', index=False)
            tag = 'AER' if aerobic else 'ANA'
            print(f'  {cond_id:24s} {tag}: bio={biomass!s:>10s}  active={active:>4d}')
            rows.append({
                'condition_id': cond_id, 'label': label, 'stage': stage,
                'O2': 'aerobic' if aerobic else 'anaerobic',
                'c_sources': ';'.join(f'{k}:{v}' for k, v in c_sources.items()),
                'status': status,
                'biomass': round(biomass, 6) if biomass is not None else None,
                'active_fluxes': active,
                'C_uptake_theor_mmolC': round(c_in_theor, 3),
                'C_uptake_realized_mmolC': round(c_in_real, 3),
                'biomass_per_C': round(bio_per_C, 6),
                'notes': notes,
            })

    pd.DataFrame(rows).to_csv(summary, sep='\t', index=False)
    print(f'\nsaved: {summary}')


if __name__ == '__main__':
    model_path = sys.argv[1] if len(sys.argv) > 1 else 'models/fsp237_gapfilled_Version1_noGenes.json'
    suffix     = sys.argv[2] if len(sys.argv) > 2 else 'v1'
    if not os.path.isabs(model_path):
        model_path = os.path.join(HERE, model_path)
    run_panel(model_path, suffix)
