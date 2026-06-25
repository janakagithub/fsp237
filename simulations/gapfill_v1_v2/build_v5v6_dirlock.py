#!/usr/bin/env python3
"""V5 / V6 -- lock fatty-acid beta-oxidation and the fungal Ashwell
D-galacturonate pathway to their degradation direction only.

Inputs:
  V3 : models/fsp237_gapfilled_Version3_dedup_noGenes.json
  V4 : models/fsp237_gapfilled_Version4_dedup_genes_integrated.json

Outputs:
  V5 : models/fsp237_gapfilled_Version5_dirlock_noGenes.json
  V6 : models/fsp237_gapfilled_Version6_dirlock_genes_integrated.json
  reports/v5_direction_locks.tsv  -- per-rxn before/after bounds + which
                                      direction in the ModelSEED-written eq
                                      corresponds to degradation.

Logic: for each enzyme in the beta-oxidation cycle and the Ashwell pathway,
we determine which sign of flux corresponds to the DEGRADATION direction
(based on the ModelSEED-written equation orientation), then set the bounds
to permit only that direction:
  - 'forward'  -> bounds = (0, ub_current_or_1000)
  - 'reverse'  -> bounds = (lb_current_or_-1000, 0)

Reactions kept reversible:
  - Peroxisomal cofactor shuttles (tx_*) -- transport stubs, biology agnostic
  - Penttilae L-arabinose pathway        -- user explicitly out of scope
  - Acyl-CoA ligases (rxn00947, rxn09445, rxn05736) -- activation can run
    either way for different substrates; not part of the chain-shortening
    cycle itself

This prevents the optimizer from using these reactions as a side-route for
biosynthesis (e.g., synthesizing palmitoyl-CoA from acetyl-CoA via reverse
beta-ox, or producing pectin-derived galacturonate from glycerol via
reverse Ashwell).
"""
import csv
import os

import cobra

BASE = '/home/janakae/fungalTemplate/imm904CobraModel'
HERE = f'{BASE}/simulations/gapfill_v1_v2'
V3   = f'{HERE}/models/fsp237_gapfilled_Version3_dedup_noGenes.json'
V4   = f'{HERE}/models/fsp237_gapfilled_Version4_dedup_genes_integrated.json'
V5   = f'{HERE}/models/fsp237_gapfilled_Version5_dirlock_noGenes.json'
V6   = f'{HERE}/models/fsp237_gapfilled_Version6_dirlock_genes_integrated.json'
LOG  = f'{HERE}/reports/v5_direction_locks.tsv'


# Per-rxn direction spec: which way of the ModelSEED equation = degradation
# 'forward' = positive flux is degradation -> bounds become (0, ub)
# 'reverse' = negative flux is degradation -> bounds become (lb, 0)
#
# Derived by reading each equation carefully (see PATHWAY_DIAGRAMS.md). Comments
# show the reaction as ModelSEED writes it; the spec column says which sign
# of flux corresponds to actual chain-shortening degradation.
DIRECTION_LOCKS = {
    # ----- Beta-oxidation: STEP 1 (acyl-CoA oxidase / dehydrogenase) -----
    # Acyl-CoA -> trans-2-enoyl-CoA; written degradation-forward.
    'rxn09474_x0': ('forward', 'C26 oxidase: acyl-CoA + O2 -> enoyl-CoA + H2O2'),
    'rxn09475_x0': ('forward', 'C18 oxidase'),
    'rxn08053_x0': ('forward', 'C18 FAD-DH: acyl-CoA + FAD -> enoyl-CoA + FADH2'),
    'rxn09476_x0': ('forward', 'C16 oxidase'),
    'rxn09477_x0': ('forward', 'C14 oxidase'),
    'rxn09478_x0': ('forward', 'C12 oxidase'),
    'rxn02720_x0': ('forward', 'C12 FAD-DH'),
    'rxn09479_x0': ('forward', 'C10 oxidase'),
    'rxn02679_x0': ('forward', 'C8 FAD-DH (V1 added)'),
    'rxn03251_x0': ('forward', 'C6 FAD-DH (V1 added)'),
    'rxn00868_x0': ('forward', 'C4 NAD-DH (V1 added): butyryl-CoA -> crotonyl-CoA'),

    # ----- Beta-oxidation: STEP 2 (enoyl-CoA hydratase) -----
    # ModelSEED writes: 3-OH-acyl-CoA <=> H2O + 2-trans-enoyl-CoA  (synthesis/
    # dehydration direction). DEGRADATION = REVERSE (enoyl + H2O -> 3-OH).
    'rxn09473_x0': ('forward', 'C26 hydratase: H2O + enoyl <-> 3-OH; written this way -> degradation FORWARD'),
    'rxn08413_x0': ('reverse', 'C18 hydratase: written as 3-OH -> H2O + enoyl, degradation reverse'),
    'rxn03240_x0': ('reverse', 'C16 hydratase'),
    'rxn03241_x0': ('reverse', 'C14 hydratase'),
    'rxn02911_x0': ('reverse', 'C12 hydratase'),
    'rxn03245_x0': ('reverse', 'C10 hydratase'),
    'rxn03247_x0': ('reverse', 'C8 hydratase (V1 added)'),
    'rxn03250_x0': ('reverse', 'C6 hydratase (V1 added)'),
    'rxn02167_x0': ('reverse', 'C4 hydratase (V1 added)'),

    # ----- Beta-oxidation: STEP 3 (3-OH-acyl-CoA dehydrogenase) -----
    # Goal: 3-OH-acyl-CoA + NAD -> 3-oxoacyl-CoA + NADH (degradation).
    # ModelSEED notation varies; we pick the right sign per reaction.
    'rxn09461_x0': ('forward', 'C26 3-OH-DH: NAD + 3-OH <-> NADH + H+ + 3-oxo; degradation forward'),
    'rxn20476_x0': ('reverse', 'C18 3-OH-DH: NADH + H+ + 3-oxo <-> NAD + 3-OH; degradation reverse'),
    'rxn20474_x0': ('reverse', 'C16 3-OH-DH'),
    'rxn30390_x0': ('reverse', 'C14 3-OH-DH'),
    'rxn20478_x0': ('reverse', 'C12 3-OH-DH'),
    'rxn20480_x0': ('reverse', 'C10 3-OH-DH'),
    'rxn03246_x0': ('forward', 'C8 3-OH-DH (V1 added): NAD + 3-OH <-> NADH + H+ + 3-oxo; forward'),
    'rxn03249_x0': ('forward', 'C6 3-OH-DH (V1 added): forward'),
    'rxn03861_x0': ('forward', 'C4 acetoacetyl-CoA reductase (V1 added, NADP)'),

    # ----- Beta-oxidation: STEP 4 (3-ketoacyl-CoA thiolase, chain shortening) -----
    # Goal: 3-oxoacyl-CoA + CoA -> acetyl-CoA + acyl-CoA(n-2).
    # rxn08767_x0 (C18->C16) is already written and bounded in the
    # degradation direction (no change needed).
    'rxn02804_x0': ('reverse', 'C16->C14 thiolase: AcCoA + C14 <-> CoA + 3-oxoC16; degradation reverse'),
    'rxn06510_x0': ('reverse', 'C14->C12 thiolase'),
    'rxn03243_x0': ('reverse', 'C12->C10 thiolase'),
    'rxn02680_x0': ('reverse', 'C10->C8 thiolase'),
    'rxn03248_x0': ('reverse', 'C8->C6 thiolase'),
    'rxn00874_x0': ('reverse', 'C6->C4 thiolase'),
    'rxn00178_x0': ('reverse', 'C4 terminal thiolase: 2 AcCoA <-> CoA + acetoacetyl-CoA; degradation reverse'),

    # ----- Ashwell D-galacturonate pathway (cytosolic) -----
    # User says fungi don't synthesize these intermediates -- lock to degradation.
    'rxn05673_c0': ('forward', 'galU/H+ symport: H+_e0 + galU_e0 <-> H+_c0 + galU_c0; only import'),
    'rxn07491_c0': ('reverse', 'GAR1: NADP + L-galactonate <-> NADPH + H+ + galU; degradation reverse (galU + NADPH -> L-galactonate)'),
    # rxn21749_c0 LGD1 already forward-only (0, 1000); no change needed.
    'rxn21750_c0': ('forward', 'LGA1: KDG <-> pyruvate + L-glyceraldehyde; degradation forward (cleavage)'),
    'rxn09954_c0': ('forward', 'GLD1: NADPH + H+ + L-glyceraldehyde <-> NADP + glycerol; degradation forward (reduce)'),
}


def apply_locks(in_path, out_path):
    print(f'\n========== {os.path.basename(in_path)} ==========')
    m = cobra.io.load_json_model(in_path)
    print(f'  start: {len(m.reactions)} rxns, {len(m.genes)} genes')

    log_rows = []
    n_changed = 0
    for rid, (direction, note) in DIRECTION_LOCKS.items():
        if rid not in [r.id for r in m.reactions]:
            log_rows.append({'rxn_id': rid, 'status': 'missing',
                              'old_bounds': '', 'new_bounds': '',
                              'direction_locked': direction, 'note': note})
            continue
        r = m.reactions.get_by_id(rid)
        old = (r.lower_bound, r.upper_bound)
        if direction == 'forward':
            # Permit only positive flux
            r.lower_bound = 0.0
            r.upper_bound = old[1] if old[1] > 0 else 1000.0
        elif direction == 'reverse':
            # Permit only negative flux
            r.lower_bound = old[0] if old[0] < 0 else -1000.0
            r.upper_bound = 0.0
        else:
            raise ValueError(f'Unknown direction "{direction}" for {rid}')
        new = (r.lower_bound, r.upper_bound)
        changed = old != new
        if changed: n_changed += 1
        log_rows.append({'rxn_id': rid, 'status': 'locked' if changed else 'already',
                          'old_bounds': str(old), 'new_bounds': str(new),
                          'direction_locked': direction, 'note': note})

    print(f'  {n_changed} reactions had bounds tightened to degradation-only')
    print(f'  end  : {len(m.reactions)} rxns, {len(m.genes)} genes')
    cobra.io.save_json_model(m, out_path)
    print(f'  saved: {out_path}')

    return log_rows


def main():
    rows_v5 = apply_locks(V3, V5)
    rows_v6 = apply_locks(V4, V6)
    # Write one combined log (the spec is the same; only bounds differ if any)
    with open(LOG, 'w', newline='') as fh:
        cols = ['rxn_id','status','old_bounds','new_bounds','direction_locked','note']
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        w.writeheader()
        for row in rows_v5: w.writerow(row)
    print(f'\nlog: {LOG}')


if __name__ == '__main__':
    main()
