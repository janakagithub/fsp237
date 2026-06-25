#!/usr/bin/env python3
"""V9 / V10 -- two cleanups applied to V5 / V6:

(1) Complete VLCFA degradation chain with biosynthesis-consistent cpd IDs.

The model's source had two broken reactions referencing an empty
metabolite `cpd15297_x0` (no name, no formula):
  rxn09480_x0  C26->C24 thiolase  (cpd15297 was the C24-CoA product)
  rxn09463_x0  lumped 3-round C24->C18 beta-oxidation (cpd15297 was substrate)

Meanwhile the biosynthesis pathway in _c0 / _r0 uses the canonical IDs:
  C24-CoA  = cpd29705   (tetracosanoyl-CoA, formula C45H78N7O17P3S)
  C22-CoA  = cpd16343   (Behenoyl-CoA)
  C20-CoA  = cpd01393   (Icosanoyl-CoA)
  + 3-OH, 3-oxo, trans-2-enoyl variants per chain

We FIX the two broken reactions by swapping cpd15297_x0 -> cpd29705_x0
(removes the orphan), then ADD the chain-specific C24/C22/C20 cycle (12
new reactions), all direction-locked to degradation. New _x0 metabolites
copy names/formulas/charges from the _r0 biosynthesis-side equivalents.

(2) Drop reactions in compartments g0 / n0 / v0.

The 88 reactions in Golgi / nucleus / vacuole compartments carry zero
flux across all 36 panel simulations -- the biomass has no organelle-
specific demand for those compartments. They're inherited yeast-specific
trafficking reactions that bloat the model. Removing them simplifies the
model without changing any biology that the simulation panel tests.

Inputs:
  V5 : models/fsp237_gapfilled_Version5_dirlock_noGenes.json
  V6 : models/fsp237_gapfilled_Version6_dirlock_genes_integrated.json

Outputs:
  V9  : models/fsp237_gapfilled_Version9_vlcfa_complete_noGenes.json
  V10 : models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json
  reports/v9_change_log.tsv
"""
import csv
import os
import re

import cobra

BASE = '/home/janakae/fungalTemplate/imm904CobraModel'
HERE = f'{BASE}/simulations/gapfill_v1_v2'
V5  = f'{HERE}/models/fsp237_gapfilled_Version5_dirlock_noGenes.json'
V6  = f'{HERE}/models/fsp237_gapfilled_Version6_dirlock_genes_integrated.json'
V9  = f'{HERE}/models/fsp237_gapfilled_Version9_vlcfa_complete_noGenes.json'
V10 = f'{HERE}/models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json'
LOG = f'{HERE}/reports/v9_change_log.tsv'

# Compartments to drop completely (no flux in any panel simulation;
# yeast-specific trafficking inherited from iMM904 source).
DROP_COMPARTMENTS = {'g0', 'n0', 'v0'}

# Metabolites we'll create in _x0 by copying names/formulas/charges from
# their biosynthesis-side _r0 equivalents. Empty entries mean already in _x0.
NEW_X0_METS = {
    # base : (preferred name, formula, charge, kegg_id)
    'cpd29705': ('tetracosanoyl-CoA (n-C24:0CoA)',    'C45H78N7O17P3S', -4, ''),  # already present
    'cpd24272': ('(2E)-lignocerenoyl-CoA',            'C45H76N7O17P3S', -4, ''),  # C24 trans-2-enoyl
    'cpd24268': ('(3R)-3-hydroxy-lignoceroyl-CoA',    'C45H78N7O18P3S', -4, ''),  # C24 3-OH
    'cpd24264': ('3-oxo-lignoceroyl-CoA',             'C45H76N7O18P3S', -4, ''),  # C24 3-oxo
    'cpd16343': ('Behenoyl-CoA',                      'C43H74N7O17P3S', -4, 'C04425'),  # C22-CoA
    'cpd24271': ('(2E)-docosenoyl-CoA',               'C43H72N7O17P3S', -4, ''),  # C22 trans-2-enoyl
    'cpd24267': ('(3R)-3-hydroxy-behenoyl-CoA',       'C43H74N7O18P3S', -4, ''),  # C22 3-OH
    'cpd22630': ('3-oxo-behenoyl-CoA',                'C43H72N7O18P3S', -4, ''),  # C22 3-oxo
    'cpd01393': ('Icosanoyl-CoA',                     'C41H70N7O17P3S', -4, 'C02232'),  # C20-CoA
    'cpd24270': ('(2E)-arachidoenoyl-CoA',            'C41H68N7O17P3S', -4, ''),  # C20 trans-2-enoyl
    'cpd24266': ('(3R)-3-hydroxy-arachidoyl-CoA',     'C41H70N7O18P3S', -4, ''),  # C20 3-OH
    'cpd24263': ('3-oxo-arachidoyl-CoA',              'C41H68N7O18P3S', -4, ''),  # C20 3-oxo
}

# 12 new VLCFA chain reactions (C24/C22/C20 cycle, degradation direction only).
# Plus 1 missing transport (hexacosanoate e0 -> c0 -- analog of palmitate's
# HDCAt_c0; without this C26 can't enter the cytosol let alone be activated).
# All four cycle steps per chain + chain-shortening thiolase to next stage.
# rxn09480_x0 (C26->C24 thiolase) is fixed in-place, not re-added.
VLCFA_RXNS = [
    # ---- Membrane transport ----
    ('HXCAt', 'Hexacosanoate transport via facilitated diffusion (e0 -> c0)',
     'cpd15240_e0 => cpd15240_c0',
     ''),
    # ---- C24 cycle ----
    ('vlcfa_oxC24', 'C24 acyl-CoA dehydrogenase (FAD)',
     'cpd00015_x0 + cpd29705_x0 => cpd00982_x0 + cpd24272_x0',
     '1.3.3.6'),
    ('vlcfa_hyC24', 'C24 enoyl-CoA hydratase',
     'cpd00001_x0 + cpd24272_x0 => cpd24268_x0',
     '4.2.1.17'),
    ('vlcfa_dhC24', 'C24 3-OH-acyl-CoA dehydrogenase (NAD)',
     'cpd00003_x0 + cpd24268_x0 => cpd00004_x0 + cpd00067_x0 + cpd24264_x0',
     '1.1.1.35'),
    ('vlcfa_thC24', '3-ketoacyl-CoA thiolase (C24->C22)',
     'cpd00010_x0 + cpd24264_x0 => cpd00022_x0 + cpd16343_x0',
     '2.3.1.16'),
    # ---- C22 cycle ----
    ('vlcfa_oxC22', 'C22 acyl-CoA dehydrogenase (FAD)',
     'cpd00015_x0 + cpd16343_x0 => cpd00982_x0 + cpd24271_x0',
     '1.3.3.6'),
    ('vlcfa_hyC22', 'C22 enoyl-CoA hydratase',
     'cpd00001_x0 + cpd24271_x0 => cpd24267_x0',
     '4.2.1.17'),
    ('vlcfa_dhC22', 'C22 3-OH-acyl-CoA dehydrogenase (NAD)',
     'cpd00003_x0 + cpd24267_x0 => cpd00004_x0 + cpd00067_x0 + cpd22630_x0',
     '1.1.1.35'),
    ('vlcfa_thC22', '3-ketoacyl-CoA thiolase (C22->C20)',
     'cpd00010_x0 + cpd22630_x0 => cpd00022_x0 + cpd01393_x0',
     '2.3.1.16'),
    # ---- C20 cycle ----
    ('vlcfa_oxC20', 'C20 acyl-CoA dehydrogenase (FAD)',
     'cpd00015_x0 + cpd01393_x0 => cpd00982_x0 + cpd24270_x0',
     '1.3.3.6'),
    ('vlcfa_hyC20', 'C20 enoyl-CoA hydratase',
     'cpd00001_x0 + cpd24270_x0 => cpd24266_x0',
     '4.2.1.17'),
    ('vlcfa_dhC20', 'C20 3-OH-acyl-CoA dehydrogenase (NAD)',
     'cpd00003_x0 + cpd24266_x0 => cpd00004_x0 + cpd00067_x0 + cpd24263_x0',
     '1.1.1.35'),
    # C20 -> C18 (stearoyl-CoA = cpd00327, alias 'strcoa', already in _x0)
    ('vlcfa_thC20', '3-ketoacyl-CoA thiolase (C20->C18)',
     'cpd00010_x0 + cpd24263_x0 => cpd00022_x0 + cpd00327_x0',
     '2.3.1.16'),
]

TERM_RE  = re.compile(r'(?:\((\d+(?:\.\d+)?)\)\s*)?(cpd\d+_[a-z]\d+)')
ARROW_RE = re.compile(r'\s*(<=>|=>|<=|->|-->|<-|<--)\s*')


def parse_eq(eq):
    arrow = ARROW_RE.search(eq).group(1)
    lhs, rhs = eq.split(arrow, 1)
    stoich = {}
    for side, sign in [(lhs, -1), (rhs, +1)]:
        for coef_str, mid in TERM_RE.findall(side):
            coef = float(coef_str) if coef_str else 1.0
            stoich[mid] = stoich.get(mid, 0.0) + sign * coef
    return stoich


def ensure_met(model, mid):
    if mid in [x.id for x in model.metabolites]:
        return model.metabolites.get_by_id(mid)
    base, comp = mid.rsplit('_', 1)
    spec = NEW_X0_METS.get(base)
    if spec:
        name, formula, charge, kegg = spec
        annot = {'modelseed.compound': base}
        if kegg: annot['kegg.compound'] = kegg
    else:
        name, formula, charge, annot = base, '', 0, {}
    met = cobra.Metabolite(id=mid, name=f'{name}_{comp}', formula=formula,
                            charge=charge, compartment=comp)
    met.annotation = annot
    model.add_metabolites([met])
    return met


def fix_broken_rxn(model, rid, swap_from, swap_to):
    """Replace `swap_from` metabolite reference with `swap_to` in reaction rid."""
    if rid not in [r.id for r in model.reactions]: return None
    r = model.reactions.get_by_id(rid)
    if swap_from not in [m.id for m in r.metabolites]: return None
    from_met = model.metabolites.get_by_id(swap_from)
    to_met = ensure_met(model, swap_to)
    coef = r.metabolites[from_met]
    r.add_metabolites({from_met: -coef, to_met: coef})
    return (rid, swap_from, swap_to, coef)


def build(in_path, out_path, label):
    print(f'\n========== {label} from {os.path.basename(in_path)} ==========')
    m = cobra.io.load_json_model(in_path)
    print(f'  start: {len(m.reactions)} rxns, {len(m.metabolites)} mets, {len(m.genes)} genes')
    log_rows = []

    # --- A. Drop reactions in unused compartments ---
    to_drop = []
    for r in m.reactions:
        comps = {met.compartment for met in r.metabolites}
        if comps & DROP_COMPARTMENTS:
            to_drop.append(r.id)
    m.remove_reactions(to_drop, remove_orphans=True)
    print(f'  dropped {len(to_drop)} reactions touching {sorted(DROP_COMPARTMENTS)}')
    for rid in to_drop:
        log_rows.append({'op': 'drop_compartment', 'rxn_id': rid, 'detail': 'g0/n0/v0 reaction unused in panel'})

    # --- B. Fix broken cpd15297_x0 references ---
    for rid in ['rxn09480_x0', 'rxn09463_x0']:
        result = fix_broken_rxn(m, rid, 'cpd15297_x0', 'cpd29705_x0')
        if result:
            print(f'  fixed {result[0]}: cpd15297_x0 -> cpd29705_x0 (coef {result[3]})')
            log_rows.append({'op': 'fix_orphan_cpd', 'rxn_id': result[0],
                              'detail': 'cpd15297_x0 (empty) -> cpd29705_x0 (canonical C24-CoA)'})

    # Remove orphan cpd15297_x0 if now unused
    if 'cpd15297_x0' in [x.id for x in m.metabolites]:
        met = m.metabolites.get_by_id('cpd15297_x0')
        if len(met.reactions) == 0:
            m.remove_metabolites([met])
            print('  removed orphan cpd15297_x0')
            log_rows.append({'op': 'drop_orphan_met', 'rxn_id': 'cpd15297_x0',
                              'detail': 'no remaining producers/consumers after fix'})

    # --- C. Add VLCFA transport + chain-specific reactions ---
    added = 0
    for base, name, eq, ec in VLCFA_RXNS:
        # transport reaction id stays bare (HXCAt); chain reactions get _x0 suffix
        rid = base if base.startswith('HXCAt') else f'{base}_x0'
        if rid in [r.id for r in m.reactions]:
            continue
        stoich = parse_eq(eq)
        met_d = {ensure_met(m, mid): coef for mid, coef in stoich.items()}
        # All VLCFA reactions are degradation-direction-only (forward)
        r = cobra.Reaction(rid, name=name, lower_bound=0.0, upper_bound=1000.0)
        r.add_metabolites(met_d)
        r.subsystem = 'VLCFA transport' if rid == 'HXCAt' else 'beta-oxidation (VLCFA)'
        r.annotation = {'ec-code': ec} if ec else {}
        m.add_reactions([r])
        added += 1
        log_rows.append({'op': 'add_rxn', 'rxn_id': rid, 'detail': f'{name}; EC {ec}'})
    print(f'  added {added} VLCFA reactions (1 transport + C24/C22/C20 cycle)')
    print(f'  end  : {len(m.reactions)} rxns, {len(m.metabolites)} mets, {len(m.genes)} genes')
    cobra.io.save_json_model(m, out_path)
    print(f'  saved: {out_path}')
    return log_rows


def main():
    rows9  = build(V5,  V9,  'V9  (no genes)')
    rows10 = build(V6,  V10, 'V10 (with genes)')
    # Log just V9 -- same ops apply for V10
    with open(LOG, 'w', newline='') as fh:
        cols = ['op', 'rxn_id', 'detail']
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        w.writeheader()
        for row in rows9: w.writerow(row)
    print(f'\nlog: {LOG}')


if __name__ == '__main__':
    main()
