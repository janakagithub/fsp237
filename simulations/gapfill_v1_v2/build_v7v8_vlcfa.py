#!/usr/bin/env python3
"""V7 / V8 -- close the VLCFA gap (hexacosanoate aerobic = 0 in V5/V6).

Adds the peroxisomal β-oxidation cycle for chain lengths C24, C22, C20
so that hexacosanoate (C26) can be fully degraded:

  C26 -> C24 -> C22 -> C20 -> C18 (existing chain handles C18 down to C4)

For each new chain length N ∈ {24, 22, 20} the standard 4-step cycle
is added (oxidase, hydratase, 3-OH-DH, chain-shortening thiolase), plus
the missing C26→C24 thiolase. All 13 new reactions are forward-only
(degradation direction) to match the V5/V6 dirlock policy.

Inputs:
  V5 : models/fsp237_gapfilled_Version5_dirlock_noGenes.json
  V6 : models/fsp237_gapfilled_Version6_dirlock_genes_integrated.json

Outputs:
  V7 : models/fsp237_gapfilled_Version7_vlcfa_noGenes.json
  V8 : models/fsp237_gapfilled_Version8_vlcfa_genes_integrated.json
  reports/v7_vlcfa_added.tsv  -- per-rxn audit
"""
import csv
import os
import re

import cobra

BASE = '/home/janakae/fungalTemplate/imm904CobraModel'
HERE = f'{BASE}/simulations/gapfill_v1_v2'
V5 = f'{HERE}/models/fsp237_gapfilled_Version5_dirlock_noGenes.json'
V6 = f'{HERE}/models/fsp237_gapfilled_Version6_dirlock_genes_integrated.json'
V7 = f'{HERE}/models/fsp237_gapfilled_Version7_vlcfa_noGenes.json'
V8 = f'{HERE}/models/fsp237_gapfilled_Version8_vlcfa_genes_integrated.json'
LOG = f'{HERE}/reports/v7_vlcfa_added.tsv'

# New metabolites needed for the VLCFA chain (ModelSEED cpd IDs where they
# exist; biology-correct names + formulas inferred from chain length).
NEW_METS = {
    # Acyl-CoA family: -2 charge for free CoA-ester carboxylate not relevant
    # here (all are thioesters with -4 net at physiological pH). Formulas
    # follow CnH(2n-3)N7O17P3S for saturated acyl-CoA, -2H for trans-2-enoyl,
    # +O for 3-OH, no extra mass for 3-oxo (loses 2H, gains O = +14).
    # We trust ModelSEED naming where available.
    # ---- C24 (Lignoceric / Tetracosanoic) ----
    # cpd15297 (Lignoceroyl-CoA) already in V5/V6
    'cpd16840': ('trans-2-Tetracosenoyl-CoA',     'C45H78N7O17P3S', -4, ''),
    'cpd16774': ('(S)-3-Hydroxytetracosanoyl-CoA','C45H80N7O18P3S', -4, ''),
    'cpd16786': ('3-Oxotetracosanoyl-CoA',        'C45H78N7O18P3S', -4, ''),
    # ---- C22 (Behenic) ----
    'cpd16343': ('Behenoyl-CoA',                  'C43H76N7O17P3S', -4, 'C04425'),
    'cpd16839': ('trans-2-Docosenoyl-CoA',        'C43H74N7O17P3S', -4, ''),
    'cpd16773': ('(S)-3-Hydroxydocosanoyl-CoA',   'C43H76N7O18P3S', -4, ''),
    'cpd16785': ('3-Oxodocosanoyl-CoA',           'C43H74N7O18P3S', -4, ''),
    # ---- C20 (Icosanoic / Arachidic) ----
    'cpd01393': ('Icosanoyl-CoA',                 'C41H72N7O17P3S', -4, 'C02232'),
    'cpd16838': ('trans-2-Icosenoyl-CoA',         'C41H70N7O17P3S', -4, ''),
    'cpd16772': ('(S)-3-Hydroxyicosanoyl-CoA',    'C41H72N7O18P3S', -4, ''),
    'cpd16784': ('3-Oxoicosanoyl-CoA',            'C41H70N7O18P3S', -4, ''),
}

# (rxn_base, label, equation, direction_after_lock, ec, subsystem)
# All NEW reactions are forward-only (degradation direction) to comply with
# V5/V6 dirlock. Where an existing ModelSEED rxn ID exists (e.g., rxn09480)
# we use it; otherwise we coin a stable custom id `vlcfa_<step>_C<chain>`.
VLCFA_RXNS = [
    # C26 -> C24 thiolase
    ('rxn09480', 'x0', '3-ketoacyl-CoA thiolase (C26->C24)',
     'cpd00010_x0 + cpd15207_x0 => cpd00022_x0 + cpd15297_x0',
     '>', '2.3.1.16', 'beta-oxidation (VLCFA chain-shortening)'),
    # ---- C24 cycle ----
    ('vlcfa_oxC24', 'x0', 'C24 acyl-CoA dehydrogenase (FAD)',
     'cpd00015_x0 + cpd15297_x0 => cpd00982_x0 + cpd16840_x0',
     '>', '1.3.3.6', 'beta-oxidation (VLCFA C24)'),
    ('vlcfa_hyC24', 'x0', 'C24 enoyl-CoA hydratase',
     'cpd00001_x0 + cpd16840_x0 => cpd16774_x0',
     '>', '4.2.1.17', 'beta-oxidation (VLCFA C24)'),
    ('vlcfa_dhC24', 'x0', 'C24 3-OH-acyl-CoA dehydrogenase (NAD)',
     'cpd00003_x0 + cpd16774_x0 => cpd00004_x0 + cpd00067_x0 + cpd16786_x0',
     '>', '1.1.1.35', 'beta-oxidation (VLCFA C24)'),
    ('vlcfa_thC24', 'x0', '3-ketoacyl-CoA thiolase (C24->C22)',
     'cpd00010_x0 + cpd16786_x0 => cpd00022_x0 + cpd16343_x0',
     '>', '2.3.1.16', 'beta-oxidation (VLCFA chain-shortening)'),
    # ---- C22 cycle ----
    ('vlcfa_oxC22', 'x0', 'C22 acyl-CoA dehydrogenase (FAD)',
     'cpd00015_x0 + cpd16343_x0 => cpd00982_x0 + cpd16839_x0',
     '>', '1.3.3.6', 'beta-oxidation (VLCFA C22)'),
    ('vlcfa_hyC22', 'x0', 'C22 enoyl-CoA hydratase',
     'cpd00001_x0 + cpd16839_x0 => cpd16773_x0',
     '>', '4.2.1.17', 'beta-oxidation (VLCFA C22)'),
    ('vlcfa_dhC22', 'x0', 'C22 3-OH-acyl-CoA dehydrogenase (NAD)',
     'cpd00003_x0 + cpd16773_x0 => cpd00004_x0 + cpd00067_x0 + cpd16785_x0',
     '>', '1.1.1.35', 'beta-oxidation (VLCFA C22)'),
    ('vlcfa_thC22', 'x0', '3-ketoacyl-CoA thiolase (C22->C20)',
     'cpd00010_x0 + cpd16785_x0 => cpd00022_x0 + cpd01393_x0',
     '>', '2.3.1.16', 'beta-oxidation (VLCFA chain-shortening)'),
    # ---- C20 cycle ----
    ('vlcfa_oxC20', 'x0', 'C20 acyl-CoA dehydrogenase (FAD)',
     'cpd00015_x0 + cpd01393_x0 => cpd00982_x0 + cpd16838_x0',
     '>', '1.3.3.6', 'beta-oxidation (VLCFA C20)'),
    ('vlcfa_hyC20', 'x0', 'C20 enoyl-CoA hydratase',
     'cpd00001_x0 + cpd16838_x0 => cpd16772_x0',
     '>', '4.2.1.17', 'beta-oxidation (VLCFA C20)'),
    ('vlcfa_dhC20', 'x0', 'C20 3-OH-acyl-CoA dehydrogenase (NAD)',
     'cpd00003_x0 + cpd16772_x0 => cpd00004_x0 + cpd00067_x0 + cpd16784_x0',
     '>', '1.1.1.35', 'beta-oxidation (VLCFA C20)'),
    # C20 -> C18 (stearoyl-CoA = cpd00327, existing in model as 'strcoa_x0')
    ('vlcfa_thC20', 'x0', '3-ketoacyl-CoA thiolase (C20->C18)',
     'cpd00010_x0 + cpd16784_x0 => cpd00022_x0 + cpd00327_x0',
     '>', '2.3.1.16', 'beta-oxidation (VLCFA chain-shortening)'),
]

TERM_RE  = re.compile(r'(?:\((\d+(?:\.\d+)?)\)\s*)?(cpd\d+_[a-z]\d+)')
ARROW_RE = re.compile(r'\s*(<=>|=>|<=|->|-->|<-|<--)\s*')
DIRECTION_TO_BOUNDS = {
    '>': (0.0, 1000.0), '=>': (0.0, 1000.0),
    '<': (-1000.0, 0.0), '<=>': (-1000.0, 1000.0),
}


def parse_equation(eq_str):
    arrow = ARROW_RE.search(eq_str).group(1)
    lhs, rhs = eq_str.split(arrow, 1)
    stoich = {}
    for side, sign in [(lhs, -1), (rhs, +1)]:
        for coef_str, mid in TERM_RE.findall(side):
            coef = float(coef_str) if coef_str else 1.0
            stoich[mid] = stoich.get(mid, 0.0) + sign * coef
    return stoich


def ensure_metabolite(model, mid):
    if mid in [x.id for x in model.metabolites]:
        return model.metabolites.get_by_id(mid)
    base, comp = mid.rsplit('_', 1)
    spec = NEW_METS.get(base)
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
    print(f'    + new met: {mid:18s} ({name})')
    return met


def add_vlcfa(model):
    log = []
    for base, comp, name, eq, direction, ec, subsystem in VLCFA_RXNS:
        rid = f'{base}_{comp}'
        if rid in [r.id for r in model.reactions]:
            log.append({'rxn_id': rid, 'status': 'already_present',
                        'name': name, 'equation': eq, 'ec': ec})
            print(f'  . {rid}: already present')
            continue
        stoich = parse_equation(eq)
        met_d = {ensure_metabolite(model, mid): coef for mid, coef in stoich.items()}
        lb, ub = DIRECTION_TO_BOUNDS[direction]
        r = cobra.Reaction(rid, name=f'{name}_{comp}', lower_bound=lb, upper_bound=ub)
        r.add_metabolites(met_d)
        r.subsystem = subsystem
        annot = {'modelseed.reaction': base} if base.startswith('rxn') else {}
        if ec: annot['ec-code'] = ec
        r.annotation = annot
        model.add_reactions([r])
        log.append({'rxn_id': rid, 'status': 'added', 'name': name,
                    'equation': eq, 'ec': ec, 'direction': direction})
        print(f'  + {rid}: {name}  bounds=({lb},{ub})')
    return log


def build(in_path, out_path):
    print(f'\n========== {os.path.basename(in_path)} ==========')
    m = cobra.io.load_json_model(in_path)
    print(f'  start: {len(m.reactions)} rxns, {len(m.metabolites)} mets, {len(m.genes)} genes')
    log = add_vlcfa(m)
    print(f'  end  : {len(m.reactions)} rxns, {len(m.metabolites)} mets')
    cobra.io.save_json_model(m, out_path)
    print(f'  saved: {out_path}')
    return log


def main():
    log7 = build(V5, V7)
    log8 = build(V6, V8)
    with open(LOG, 'w', newline='') as fh:
        cols = ['rxn_id','status','name','equation','ec','direction']
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        w.writeheader()
        for row in log7:
            w.writerow({k: row.get(k, '') for k in cols})
    print(f'\nlog: {LOG}')


if __name__ == '__main__':
    main()
