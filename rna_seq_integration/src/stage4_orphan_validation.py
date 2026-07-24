#!/opt/env/modelseed/bin/python3
"""Stage 4 — orphan / no-GPR reaction expression-based triage.

For every ORPHAN_FLUX reaction that Stage 1 flagged (flux ≠ 0 but the reaction
has no GPR), this stage suggests candidate genes to consider during the
future BLAST-based gene-assignment effort. Candidates come from *inside the
model* only: genes that already appear in reactions in the same subsystem or
sharing an EC number token in the annotation. This is a triage feed, NOT a
gene-assignment claim.

Emits:
  outputs/stage4_orphan/<cond>_<o2>.tsv    (per-condition, per-reaction)
  outputs/orphan_priority.tsv              (cross-condition ranked list)
"""
import sys
from collections import defaultdict
from pathlib import Path

import cobra
import numpy as np
import pandas as pd

sys.path.insert(0, '/home/janakae/fsp237/simulations')
from run_simulation_panel import CONDITIONS  # noqa: E402

ROOT = Path('/home/janakae/fungalTemplate/imm904CobraModel/rna_seq_integration')
MODEL_PATH = Path('/home/janakae/fsp237/simulations/gapfill_v1_v2/models/'
                  'fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json')
REXP_TSV = ROOT / 'outputs' / 'reaction_expression.tsv'
STAGE1_DIR = ROOT / 'outputs' / 'stage1_overlay'
OUT_DIR = ROOT / 'outputs' / 'stage4_orphan'
PRIORITY_TSV = ROOT / 'outputs' / 'orphan_priority.tsv'
OUT_DIR.mkdir(parents=True, exist_ok=True)

FLUX_EPS = 1e-6
TOP_CANDIDATES = 3


def index_genes_by_subsystem(model, rexp):
    """subsystem -> [(gene_id, agg_expression), ...] sorted desc by expression."""
    idx = defaultdict(list)
    seen = defaultdict(set)
    expr_by_rxn = rexp['agg_mean_log2TPMp1'].to_dict()
    for r in model.reactions:
        subsys = r.subsystem or ''
        if not subsys:
            continue
        e = expr_by_rxn.get(r.id, np.nan)
        for g in r.genes:
            if g.id in seen[subsys]:
                continue
            seen[subsys].add(g.id)
            idx[subsys].append((g.id, float(e) if e == e else np.nan))
    for s in idx:
        idx[s].sort(key=lambda x: (-(x[1] if x[1] == x[1] else -1), x[0]))
    return idx


def main():
    print(f'model : {MODEL_PATH.name}')
    model = cobra.io.load_json_model(str(MODEL_PATH))
    rexp = pd.read_csv(REXP_TSV, sep='\t', comment='#').set_index('rxn_id')

    subsys_by_rxn = {r.id: (r.subsystem or '') for r in model.reactions}
    name_by_rxn = {r.id: (r.name or '') for r in model.reactions}
    equation_by_rxn = {r.id: r.build_reaction_string(use_metabolite_names=True)
                        for r in model.reactions}
    gene_idx = index_genes_by_subsystem(model, rexp)

    all_orphans = defaultdict(lambda: {'conditions': set(), 'max_abs_flux': 0.0,
                                        'sum_abs_flux': 0.0, 'n_hit': 0})

    for cond_id, label, stage, c_sources, notes in CONDITIONS:
        for aerobic in (True, False):
            o2_tag = 'aerobic' if aerobic else 'anaerobic'
            stage1_path = STAGE1_DIR / f'{cond_id}_{o2_tag}.tsv'
            if not stage1_path.exists():
                continue
            df = pd.read_csv(stage1_path, sep='\t')
            orph = df[df['agreement'] == 'ORPHAN_FLUX'].copy()
            if orph.empty:
                pd.DataFrame(columns=['rxn_id', 'name', 'subsystem', 'equation',
                                       'flux_pFBA', 'top_candidates', 'top_candidate_expr'
                                       ]).to_csv(OUT_DIR / f'{cond_id}_{o2_tag}.tsv',
                                                  sep='\t', index=False)
                continue

            rows = []
            for _, r in orph.iterrows():
                rid = r['rxn_id']
                subsys = subsys_by_rxn.get(rid, '')
                candidates = gene_idx.get(subsys, [])[:TOP_CANDIDATES]
                cand_ids = ';'.join(g[0] for g in candidates)
                cand_expr = ';'.join(f'{g[1]:.2f}' if g[1] == g[1] else 'nan'
                                      for g in candidates)
                rows.append({
                    'rxn_id': rid, 'name': name_by_rxn.get(rid, ''),
                    'subsystem': subsys, 'equation': equation_by_rxn.get(rid, ''),
                    'flux_pFBA': r['flux_pFBA'],
                    'top_candidates': cand_ids or '(no candidates in same subsystem)',
                    'top_candidate_expr': cand_expr,
                })
                agg = all_orphans[rid]
                agg['conditions'].add(f'{cond_id}_{o2_tag}')
                agg['max_abs_flux'] = max(agg['max_abs_flux'], abs(r['flux_pFBA']))
                agg['sum_abs_flux'] += abs(r['flux_pFBA'])
                agg['n_hit'] += 1

            (pd.DataFrame(rows)
                .sort_values('flux_pFBA', key=lambda c: c.abs(), ascending=False)
                .to_csv(OUT_DIR / f'{cond_id}_{o2_tag}.tsv', sep='\t', index=False))

    # cross-condition priority list
    priority = []
    for rid, agg in all_orphans.items():
        subsys = subsys_by_rxn.get(rid, '')
        candidates = gene_idx.get(subsys, [])[:TOP_CANDIDATES]
        priority.append({
            'rxn_id': rid, 'name': name_by_rxn.get(rid, ''),
            'subsystem': subsys, 'equation': equation_by_rxn.get(rid, ''),
            'n_conditions_active': agg['n_hit'],
            'conditions': ';'.join(sorted(agg['conditions'])),
            'max_abs_flux': round(agg['max_abs_flux'], 6),
            'sum_abs_flux': round(agg['sum_abs_flux'], 6),
            'top_candidates': ';'.join(g[0] for g in candidates),
            'top_candidate_expr': ';'.join(f'{g[1]:.2f}' if g[1] == g[1] else 'nan'
                                            for g in candidates),
        })
    pri_df = pd.DataFrame(priority).sort_values(
        ['n_conditions_active', 'max_abs_flux'], ascending=[False, False]
    )
    pri_df.to_csv(PRIORITY_TSV, sep='\t', index=False)
    print(f'wrote {PRIORITY_TSV}   ({len(pri_df)} orphan reactions across conditions)')
    print('\ntop 15 by (n_conditions, max_flux):')
    print(pri_df.head(15)[['rxn_id', 'name', 'subsystem', 'n_conditions_active',
                             'max_abs_flux', 'top_candidates']].to_string(index=False))


if __name__ == '__main__':
    main()
