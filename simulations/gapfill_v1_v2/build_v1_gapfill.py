#!/usr/bin/env python3
"""Phase 1 -- build FSP237 gapfilled Version 1 (no genes), all reactions in
their biologically correct compartments.

Per user instruction (2026-06-22): keep fatty-acid beta-oxidation in the
peroxisome (_x0). The existing model already has acyl-CoA oxidase + 3-OH-DH
reactions for C18, C16, C14, C12, C10, C26 chain lengths, but is missing
the shorter-chain variants and most chain-shortening thiolases, plus the
peroxisomal cofactor (ATP, NAD, AMP, PPi, CoA) shuttle reactions that the
cycle needs to actually carry flux. We add the missing chain-length-specific
reactions in _x0, plus minimal peroxisomal cofactor transport.

Gap-fill blocks added (all stoichiometry from ModelSEED):

A. Peroxisomal beta-oxidation (chain-length-specific completion)
  A1. acyl-CoA dehydrogenases (FAD-coupled, consistent with existing _x0 set):
      - rxn02679_x0  Octanoyl-CoA + FAD <-> FADH2 + (2E)-Octenoyl-CoA   (C8)
      - rxn03251_x0  Hexanoyl-CoA + FAD <-> FADH2 + (2E)-Hexenoyl-CoA   (C6)
      - rxn00868_x0  Crotonyl-CoA + NADH + H+ <-> Butyryl-CoA + NAD     (C4)
                    (only NAD-coupled variant in MS; runs reverse for beta-ox)

  A2. (2E)-Enoyl-CoA hydratases (already present for C16/C14/C12/C10):
      - rxn03247_x0  (S)-Hydroxyoctanoyl-CoA <-> H2O + (2E)-Octenoyl-CoA   (C8)
      - rxn03250_x0  (S)-Hydroxyhexanoyl-CoA <-> H2O + (2E)-Hexenoyl-CoA   (C6)
      - rxn02167_x0  (S)-3-Hydroxybutanoyl-CoA <-> H2O + Crotonyl-CoA      (C4)

  A3. (S)-3-Hydroxyacyl-CoA dehydrogenases (NAD-coupled):
      - rxn03246_x0  (S)-Hydroxyoctanoyl-CoA + NAD <-> NADH + H+ + 3-Oxooctanoyl-CoA  (C8)
      - rxn03249_x0  (S)-Hydroxyhexanoyl-CoA + NAD <-> NADH + H+ + 3-Oxohexanoyl-CoA  (C6)
      - rxn03861_x0  (S)-3-Hydroxybutanoyl-CoA + NADP <-> NADPH + H+ + Acetoacetyl-CoA (C4)

  A4. Chain-shortening 3-ketoacyl-CoA thiolases (running in beta-ox direction):
      - rxn02804_x0  (C16->C14) 3-Oxopalmitoyl-CoA + CoA -> Acetyl-CoA + Myristoyl-CoA
      - rxn06510_x0  (C14->C12) 3-Oxomyristoyl-CoA + CoA -> Acetyl-CoA + Lauroyl-CoA
      - rxn03243_x0  (C12->C10) 3-Oxolauroyl-CoA + CoA -> Acetyl-CoA + Decanoyl-CoA
      - rxn02680_x0  (C10->C8)  3-Oxodecanoyl-CoA + CoA -> Acetyl-CoA + Octanoyl-CoA
      - rxn03248_x0  (C8->C6)   3-Oxooctanoyl-CoA + CoA -> Acetyl-CoA + Hexanoyl-CoA
      - rxn00874_x0  (C6->C4)   3-Oxohexanoyl-CoA + CoA -> Acetyl-CoA + Butyryl-CoA
      - rxn00178_x0  (C4->2 acetyl-CoA) Acetoacetyl-CoA + CoA -> 2 Acetyl-CoA

  A5. Peroxisomal cofactor shuttles (NEW transport stubs):
      Without these, ATP_x0 / AMP_x0 / PPi_x0 / CoA_x0 pools are closed loops
      and cycle flux is zero. Added as reversible diffusion stubs; could be
      refined to specific antiporters (Ant1, PXA-family) in later versions.
      - tx_atp_xc      cpd00002_c0 + cpd00018_x0 <-> cpd00002_x0 + cpd00018_c0  (ATP/AMP antiport, Ant1-like)
      - tx_ppi_xc      cpd00012_x0 <-> cpd00012_c0  (PPi export)
      - tx_coa_xc      cpd00010_c0 <-> cpd00010_x0  (CoA equilibration)
      - tx_nad_xc      cpd00003_c0 <-> cpd00003_x0  (NAD/NADH peroxisomal shuttle, simplified)
      - tx_nadh_xc     cpd00004_c0 <-> cpd00004_x0
      - tx_h_xc        cpd00067_c0 <-> cpd00067_x0

B. Penttilae fungal L-arabinose pathway closure (2 cytosolic reactions)
  - rxn01391_c0  L-arabitol -> L-xylulose (EC 1.1.1.12)
  - rxn33066_c0  L-xylulose -> xylitol    (EC 1.1.1.10)

C. Fungal D-galacturonate (Ashwell) pathway (5 cytosolic reactions)
  - rxn05673_c0  D-galacturonate proton symport (e0 <-> c0)
  - rxn07491_c0  D-galU + NADPH -> L-galactonate + NADP   (GAR1, EC 1.1.1.365)
  - rxn21749_c0  L-galactonate -> 2-keto-3-deoxy-L-galactonate + H2O   (LGD1)
  - rxn21750_c0  2-keto-3-deoxy-L-galactonate -> pyruvate + L-glyceraldehyde (LGA1)
  - rxn09954_c0  L-glyceraldehyde + NADPH -> glycerol + NADP    (GLD1)

Total: 23 reactions added across compartments x0/c0.
No GPRs in this phase; gene candidates wired in Phase 4 (Version 2).

Save: simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version1_noGenes.json
"""
import json
import os
import re

import cobra

BASE = '/home/janakae/fungalTemplate/imm904CobraModel'
INPUT_MODEL  = f'{BASE}/gpr-update/fsp237_atp_safe_gsm_gpr_updated.json'
OUT_DIR      = f'{BASE}/simulations/gapfill_v1_v2'
OUT_MODEL    = f'{OUT_DIR}/models/fsp237_gapfilled_Version1_noGenes.json'
ADDED_TSV    = f'{OUT_DIR}/reports/v1_added_reactions.tsv'

os.makedirs(os.path.dirname(OUT_MODEL), exist_ok=True)
os.makedirs(os.path.dirname(ADDED_TSV), exist_ok=True)

# (rxn_id_base, compartment_suffix, name, equation, direction, ec, subsystem, reason)
GAPFILL_RXNS = [
    # ---- A1. Acyl-CoA dehydrogenases (FAD-coupled, peroxisomal) ----
    ('rxn02679', 'x0', 'Octanoyl-CoA:oxygen 2-oxidoreductase (C8)',
     'cpd00015_x0 + cpd01335_x0 <=> cpd00982_x0 + cpd03130_x0',
     '<=>', '1.3.3.6', 'beta-oxidation (acyl-CoA oxidase)',
     'C8 acyl-CoA dehydrogenase; FAD-coupled to match existing C12/C14 pattern'),
    ('rxn03251', 'x0', 'Hexanoyl-CoA:(acceptor) 2,3-oxidoreductase (C6)',
     'cpd00015_x0 + cpd03124_x0 <=> cpd00982_x0 + cpd03125_x0',
     '<=>', '1.3.3.6', 'beta-oxidation (acyl-CoA oxidase)',
     'C6 acyl-CoA dehydrogenase; FAD-coupled to match existing pattern'),
    ('rxn00868', 'x0', 'Butanoyl-CoA:(acceptor) 2,3-oxidoreductase (C4)',
     'cpd00120_x0 + cpd00003_x0 <=> cpd00650_x0 + cpd00004_x0 + cpd00067_x0',
     '<=>', '1.3.99.2', 'beta-oxidation (acyl-CoA dehydrogenase)',
     'C4 acyl-CoA dehydrogenase; NAD-coupled (only variant in MS for butyryl)'),

    # ---- A2. (2E)-Enoyl-CoA hydratases ----
    ('rxn03247', 'x0', '(S)-Hydroxyoctanoyl-CoA hydro-lyase (C8)',
     'cpd03120_x0 <=> cpd00001_x0 + cpd03130_x0',
     '<=>', '4.2.1.17', 'beta-oxidation (enoyl-CoA hydratase)',
     'C8 hydratase; was a documented gap'),
    ('rxn03250', 'x0', '(S)-Hydroxyhexanoyl-CoA hydro-lyase (C6)',
     'cpd03122_x0 <=> cpd00001_x0 + cpd03125_x0',
     '<=>', '4.2.1.17', 'beta-oxidation (enoyl-CoA hydratase)',
     'C6 hydratase; was a documented gap'),
    ('rxn02167', 'x0', '(S)-3-Hydroxybutanoyl-CoA hydro-lyase (C4)',
     'cpd00842_x0 <=> cpd00001_x0 + cpd00650_x0',
     '<=>', '4.2.1.55', 'beta-oxidation (enoyl-CoA hydratase)',
     'C4 (crotonyl-CoA) hydratase; final beta-ox cycle hydration step'),

    # ---- A3. (S)-3-Hydroxyacyl-CoA dehydrogenases ----
    ('rxn03246', 'x0', '(S)-hydroxyoctanoyl-CoA:NAD+ oxidoreductase (C8)',
     'cpd00003_x0 + cpd03120_x0 <=> cpd00004_x0 + cpd00067_x0 + cpd03121_x0',
     '<=>', '1.1.1.35', 'beta-oxidation (3-OH-acyl-CoA DH)',
     'C8 3-OH-acyl-CoA dehydrogenase'),
    ('rxn03249', 'x0', '(S)-hydroxyhexanoyl-CoA:NAD+ oxidoreductase (C6)',
     'cpd00003_x0 + cpd03122_x0 <=> cpd00004_x0 + cpd00067_x0 + cpd03123_x0',
     '<=>', '1.1.1.35', 'beta-oxidation (3-OH-acyl-CoA DH)',
     'C6 3-OH-acyl-CoA dehydrogenase'),
    ('rxn03861', 'x0', '(S)-3-Hydroxybutanoyl-CoA:NADP+ oxidoreductase (C4)',
     'cpd00006_x0 + cpd03043_x0 <=> cpd00005_x0 + cpd00067_x0 + cpd00279_x0',
     '<=>', '1.1.1.36', 'beta-oxidation (acetoacetyl-CoA reductase)',
     'Acetoacetyl-CoA reductase (NADP-coupled); C4 3-OH step'),

    # ---- A4. Chain-shortening thiolases (running in beta-ox direction = reverse of MS name) ----
    ('rxn02804', 'x0', 'myristoyl-CoA:acetyl-CoA C-myristoyltransferase (C16->C14)',
     'cpd00022_x0 + cpd01695_x0 <=> cpd00010_x0 + cpd03114_x0',
     '<=>', '2.3.1.16', 'beta-oxidation (thiolase)',
     'C16->C14 thiolase; in BO direction: 3-oxoC16 + CoA -> AcCoA + C14-CoA'),
    ('rxn06510', 'x0', 'Lauroyl-CoA:acetyl-CoA C-acyltransferase (C14->C12)',
     'cpd00022_x0 + cpd01260_x0 <=> cpd00010_x0 + cpd12689_x0',
     '<=>', '2.3.1.16', 'beta-oxidation (thiolase)',
     'C14->C12 thiolase'),
    ('rxn03243', 'x0', 'Decanoyl-CoA:acetyl-CoA C-acyltransferase (C12->C10)',
     'cpd00022_x0 + cpd03128_x0 <=> cpd00010_x0 + cpd03117_x0',
     '<=>', '2.3.1.16', 'beta-oxidation (thiolase)',
     'C12->C10 thiolase'),
    ('rxn02680', 'x0', 'Octanoyl-CoA:acetyl-CoA C-acyltransferase (C10->C8)',
     'cpd00022_x0 + cpd01335_x0 <=> cpd00010_x0 + cpd03119_x0',
     '<=>', '2.3.1.16', 'beta-oxidation (thiolase)',
     'C10->C8 thiolase'),
    ('rxn03248', 'x0', 'Hexanoyl-CoA:acetyl-CoA C-acyltransferase (C8->C6)',
     'cpd00022_x0 + cpd03124_x0 <=> cpd00010_x0 + cpd03121_x0',
     '<=>', '2.3.1.16', 'beta-oxidation (thiolase)',
     'C8->C6 thiolase'),
    ('rxn00874', 'x0', 'butanoyl-CoA:acetyl-CoA C-butanoyltransferase (C6->C4)',
     'cpd00022_x0 + cpd00120_x0 <=> cpd00010_x0 + cpd03123_x0',
     '<=>', '2.3.1.16', 'beta-oxidation (thiolase)',
     'C6->C4 thiolase'),
    ('rxn00178', 'x0', 'Acetyl-CoA:acetyl-CoA C-acetyltransferase (C4->2 AcCoA)',
     '(2) cpd00022_x0 <=> cpd00010_x0 + cpd00279_x0',
     '<=>', '2.3.1.9', 'beta-oxidation (thiolase, terminal)',
     'Final thiolase: acetoacetyl-CoA + CoA -> 2 acetyl-CoA (reversed)'),

    # ---- A5. Peroxisomal cofactor shuttles (transport stubs) ----
    ('tx_atp_xc', 'tx', 'ATP/AMP antiport (Ant1-like) peroxisomal',
     'cpd00002_c0 + cpd00018_x0 <=> cpd00002_x0 + cpd00018_c0',
     '<=>', '', 'peroxisomal cofactor transport',
     'ATP/AMP antiport stub; required for fatty acyl-CoA activation'),
    ('tx_ppi_xc', 'tx', 'PPi peroxisomal export',
     'cpd00012_x0 <=> cpd00012_c0',
     '<=>', '', 'peroxisomal cofactor transport',
     'PPi release from CoA-ligase reactions'),
    ('tx_coa_xc', 'tx', 'CoA peroxisomal equilibration',
     'cpd00010_c0 <=> cpd00010_x0',
     '<=>', '', 'peroxisomal cofactor transport',
     'Free CoA equilibration; supports activation + thiolase'),
    ('tx_nad_xc', 'tx', 'NAD peroxisomal equilibration',
     'cpd00003_c0 <=> cpd00003_x0',
     '<=>', '', 'peroxisomal cofactor transport',
     'NAD equilibration for 3-OH-acyl-CoA dehydrogenases'),
    ('tx_nadh_xc', 'tx', 'NADH peroxisomal equilibration',
     'cpd00004_c0 <=> cpd00004_x0',
     '<=>', '', 'peroxisomal cofactor transport',
     'NADH balance / shuttle equivalent'),
    ('tx_h_xc', 'tx', 'H+ peroxisomal equilibration',
     'cpd00067_c0 <=> cpd00067_x0',
     '<=>', '', 'peroxisomal cofactor transport',
     'Proton balance'),

    # ---- B. Penttilae L-arabinose pathway ----
    ('rxn01391', 'c0', 'L-arabitol:NAD+ 4-oxidoreductase (L-xylulose-forming)',
     'cpd00003_c0 + cpd00417_c0 <=> cpd00004_c0 + cpd00067_c0 + cpd00261_c0',
     '<=>', '1.1.1.12', 'L-arabinose assimilation (Penttilae)',
     'Penttilae step 2 (L-arabitol -> L-xylulose)'),
    ('rxn33066', 'c0', 'L-xylulose reductase',
     'cpd00006_c0 + cpd00306_c0 <=> cpd00005_c0 + cpd00067_c0 + cpd00261_c0',
     '<=>', '1.1.1.10', 'L-arabinose assimilation (Penttilae)',
     'Penttilae step 3 (L-xylulose -> xylitol, reversible)'),

    # ---- C. Fungal D-galacturonate (Ashwell) pathway ----
    ('rxn05673', 'c0', 'D-galacturonate transport via proton symport',
     'cpd00067_e0 + cpd00280_e0 <=> cpd00067_c0 + cpd00280_c0',
     '<=>', '', 'galU transport',
     'galU was unable to enter the cytosol; needed before any catabolism'),
    ('rxn07491', 'c0', 'D-galacturonate reductase (NADPH; GAR1, fungal)',
     'cpd00006_c0 + cpd14659_c0 <=> cpd00005_c0 + cpd00067_c0 + cpd00280_c0',
     '<=>', '1.1.1.365', 'fungal Ashwell pathway',
     'Ashwell entry (galU -> L-galactonate); NADPH-coupled fungal GAR1'),
    ('rxn21749', 'c0', 'L-galactonate dehydratase (LGD1)',
     'cpd14659_c0 => cpd00001_c0 + cpd23364_c0',
     '>', '4.2.1.146', 'fungal Ashwell pathway',
     'Ashwell step 2 (L-galactonate -> 2-keto-3-deoxy-L-galactonate)'),
    ('rxn21750', 'c0', '2-keto-3-deoxy-L-galactonate aldolase (LGA1)',
     'cpd23364_c0 <=> cpd00020_c0 + cpd01605_c0',
     '<=>', '4.1.2.55', 'fungal Ashwell pathway',
     'Ashwell step 3 (KDG -> pyruvate + L-glyceraldehyde)'),
    ('rxn09954', 'c0', 'L-glyceraldehyde reductase (GLD1)',
     'cpd00005_c0 + cpd00067_c0 + cpd01605_c0 <=> cpd00006_c0 + cpd00100_c0',
     '<=>', '1.1.1.21', 'fungal Ashwell pathway',
     'Ashwell step 4 (L-glyceraldehyde + NADPH -> glycerol)'),
]

NEW_METABOLITE_SPEC = {
    # beta-ox cycle intermediates
    'cpd00120': ('Butyryl-CoA',                   'C25H38N7O17P3S', -4, 'C00136'),
    'cpd00279': ('Acetoacetyl-CoA',               'C25H36N7O18P3S', -4, 'C00332'),
    'cpd00650': ('Crotonyl-CoA',                  'C25H38N7O17P3S', -4, 'C00877'),
    'cpd00842': ('(S)-3-Hydroxybutanoyl-CoA',     'C25H40N7O18P3S', -4, 'C01144'),
    'cpd03043': ('(S)-3-Hydroxybutyryl-CoA',      'C25H40N7O18P3S', -4, 'C01144'),
    'cpd03114': ('3-Oxopalmitoyl-CoA',            'C37H60N7O18P3S', -4, 'C05258'),
    'cpd03117': ('3-Oxododecanoyl-CoA',           'C33H52N7O18P3S', -4, 'C05262'),
    'cpd03119': ('3-Oxodecanoyl-CoA',             'C31H48N7O18P3S', -4, 'C05264'),
    'cpd03120': ('(S)-Hydroxyoctanoyl-CoA',       'C29H46N7O18P3S', -4, 'C05266'),
    'cpd03121': ('3-Oxooctanoyl-CoA',             'C29H44N7O18P3S', -4, 'C05267'),
    'cpd03122': ('(S)-Hydroxyhexanoyl-CoA',       'C27H42N7O18P3S', -4, 'C05269'),
    'cpd03123': ('3-Oxohexanoyl-CoA',             'C27H40N7O18P3S', -4, 'C05270'),
    'cpd03124': ('Hexanoyl-CoA',                  'C27H42N7O17P3S', -4, 'C05270'),
    'cpd03125': ('(2E)-Hexenoyl-CoA',             'C27H40N7O17P3S', -4, 'C05272'),
    'cpd03130': ('(2E)-Octenoyl-CoA',             'C29H44N7O17P3S', -4, 'C05276'),
    'cpd01335': ('Octanoyl-CoA',                  'C29H46N7O17P3S', -4, 'C01944'),
    'cpd01260': ('Lauroyl-CoA',                   'C33H54N7O17P3S', -4, 'C01832'),
    'cpd01695': ('Myristoyl-CoA',                 'C35H58N7O17P3S', -4, 'C02593'),
    'cpd03128': ('Decanoyl-CoA',                  'C31H50N7O17P3S', -4, 'C05274'),
    'cpd12689': ('3-Oxotetradecanoyl-CoA',        'C35H56N7O18P3S', -4, ''),
    # Penttilae + Ashwell
    'cpd00261': ('L-Lyxulose',                    'C5H10O5',         0, 'C00312'),
    'cpd00417': ('L-Arabitol',                    'C5H12O5',         0, 'C00532'),
    'cpd01605': ('L-Glyceraldehyde',              'C3H6O3',          0, 'C02426'),
    'cpd14659': ('L-Galactonate',                 'C6H11O7',        -1, 'C03383'),
    'cpd23364': ('2-keto-3-deoxy-L-galactonate',  'C6H9O6',         -1, 'C20566'),
}

TERM_RE  = re.compile(r'(?:\((\d+(?:\.\d+)?)\)\s*)?(cpd\d+_[a-z]\d+)')
ARROW_RE = re.compile(r'\s*(<=>|=>|<=|->|-->|<-|<--)\s*')

DIRECTION_TO_BOUNDS = {
    '>': (0.0, 1000.0), '=>': (0.0, 1000.0),
    '<': (-1000.0, 0.0), '<=': (-1000.0, 0.0),
    '<=>': (-1000.0, 1000.0), '<>': (-1000.0, 1000.0),
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
    if mid in [m.id for m in model.metabolites]:
        return model.metabolites.get_by_id(mid)
    base, comp = mid.rsplit('_', 1)
    spec = NEW_METABOLITE_SPEC.get(base)
    if spec:
        name, formula, charge, kegg = spec
        annot = {'modelseed.compound': base}
        if kegg: annot['kegg.compound'] = kegg
    else:
        name, formula, charge, annot = mid, '', 0, {}
    m = cobra.Metabolite(id=mid, name=f'{name}_{comp}' if not name.endswith(f'_{comp}') else name,
                          formula=formula, charge=charge, compartment=comp)
    m.annotation = annot
    model.add_metabolites([m])
    print(f'    created new metabolite: {mid:18s} ({name})')
    return m


def add_reaction(model, base, comp, name, eq, direction, ec, subsystem, reason):
    rid = f'{base}_{comp}' if comp != 'tx' else base
    if rid in [r.id for r in model.reactions]:
        return ('already_present', rid)
    stoich = parse_equation(eq)
    met_d = {ensure_metabolite(model, mid): coef for mid, coef in stoich.items()}
    lb, ub = DIRECTION_TO_BOUNDS[direction]
    r = cobra.Reaction(rid, name=name if comp == 'tx' else f'{name}_{comp}',
                        lower_bound=lb, upper_bound=ub)
    r.add_metabolites(met_d)
    r.subsystem = subsystem
    annot = {'modelseed.reaction': base} if not base.startswith('tx_') else {}
    if ec: annot['ec-code'] = ec
    r.annotation = annot
    # Intentionally NO GPR in Version 1
    model.add_reactions([r])
    return ('added', rid)


def main():
    print(f'loading: {INPUT_MODEL}')
    m = cobra.io.load_json_model(INPUT_MODEL)
    print(f'  start: {len(m.reactions)} rxns, {len(m.metabolites)} mets, {len(m.genes)} genes')

    print(f'\n=== adding {len(GAPFILL_RXNS)} gapfill reactions ===')
    log = []
    for base, comp, name, eq, direction, ec, subsystem, reason in GAPFILL_RXNS:
        status, rid = add_reaction(m, base, comp, name, eq, direction, ec, subsystem, reason)
        flag = '+' if status == 'added' else '.'
        print(f'  {flag} {rid:18s} [{status:15s}] {name[:60]}')
        log.append({'rxn_id': rid, 'status': status, 'name': name,
                    'compartment': comp, 'equation': eq, 'direction': direction,
                    'ec': ec, 'subsystem': subsystem, 'reason': reason,
                    'gpr': '(none - assigned in Version 2)'})

    print(f'\n  model now: {len(m.reactions)} rxns, {len(m.metabolites)} mets')

    cobra.io.save_json_model(m, OUT_MODEL)
    print(f'\nsaved: {OUT_MODEL}')

    import csv
    with open(ADDED_TSV, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(log[0].keys()), delimiter='\t')
        w.writeheader()
        for row in log: w.writerow(row)
    print(f'saved: {ADDED_TSV}')


if __name__ == '__main__':
    main()
