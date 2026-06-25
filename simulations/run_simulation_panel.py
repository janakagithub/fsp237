#!/usr/bin/env python3
"""Run the 18-condition FBA panel from INFECTION_SIM_PLAN.md.

Source model: gpr-update/fsp237_atp_safe_gsm_gpr_updated.json (GPR-updated
biomass-extended GSM).

For each condition we:
  - reset every exchange to lower=0  (closed)
  - re-open the standard inorganic minimal-media set (H2O, Pi, CO2, H+, NH3,
    sulfate, Fe2+, Na+, K+) at -1000
  - open the condition-specific C-source(s) at the rate spec'd in the report
  - set O2 according to aerobic / anaerobic
  - optimize bio_gsm  (and separately bio1 with glc=-1 for ATP yield)
  - record biomass, ATP yield, active-flux count, and a top-flux summary

Outputs (under imm904CobraModel/simulations/):
  - simulation_results.tsv     -- one row per condition x condition_O2
  - per_condition/<id>_aerobic.tsv  / <id>_anaerobic.tsv -- non-zero fluxes
"""
import os
import re
from collections import OrderedDict

import cobra
import pandas as pd

BASE = '/home/janakae/fungalTemplate/imm904CobraModel'
# Source the latest gapfilled + dedup + direction-locked + gene-integrated
# model (V6). This is the canonical "current" model for publication; previous
# versions are preserved on disk for reversibility (see
# simulations/gapfill_v1_v2/reports/SUMMARY.md).
MODEL_PATH = f'{BASE}/simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json'
OUT_DIR    = f'{BASE}/simulations'
PER_DIR    = f'{OUT_DIR}/per_condition'
SUMMARY    = f'{OUT_DIR}/simulation_results.tsv'

os.makedirs(PER_DIR, exist_ok=True)

# Standard inorganic minimal-media set (matches extend_biomass.py's
# media_atp_yield()).  O2 is added conditionally per simulation.
INORGANIC = {
    'cpd00001_e0': 'H2O',
    'cpd00009_e0': 'Phosphate',
    'cpd00011_e0': 'CO2',
    'cpd00067_e0': 'H+',
    'cpd00013_e0': 'NH3',
    'cpd00048_e0': 'Sulfate',
    'cpd10515_e0': 'Fe2+',
    'cpd00971_e0': 'Na+',
    'cpd00205_e0': 'K+',
}
O2 = 'cpd00007_e0'

# All 20 proteinogenic amino-acid exchanges in the model
AA_EXCHANGES = [
    'cpd00023_e0',  # Glu
    'cpd00033_e0',  # Gly
    'cpd00035_e0',  # Ala
    'cpd00039_e0',  # Lys
    'cpd00041_e0',  # Asp
    'cpd00051_e0',  # Arg
    'cpd00053_e0',  # Gln
    'cpd00054_e0',  # Ser
    'cpd00060_e0',  # Met
    'cpd00065_e0',  # Trp
    'cpd00066_e0',  # Phe
    'cpd00069_e0',  # Tyr
    'cpd00084_e0',  # Cys
    'cpd00107_e0',  # Leu
    'cpd00119_e0',  # His
    'cpd00129_e0',  # Pro
    'cpd00132_e0',  # Asn
    'cpd00156_e0',  # Val
    'cpd00161_e0',  # Thr
    'cpd00322_e0',  # Ile
]

# (id, label, stage, C-source dict {met_id: -uptake}, notes)
# Negative numbers are uptake rates (mmol/gDW/h); positive in the exchange
# convention means production. We pass the magnitude only and apply the
# sign when setting lower_bound.
CONDITIONS = [
    # ---- Pre-infection ----
    ('01_trehalose',     'Trehalose',          'pre-infection', {'cpd00794_e0': 1.0}, 'Spore-stored TRHL mobilization'),
    ('02_glycerol',      'Glycerol',           'pre-infection', {'cpd00100_e0': 1.0}, 'Appressorial osmolyte / lipid backbone'),
    ('03_palmitate',     'Palmitate (C16:0)',  'pre-infection', {'cpd00214_e0': 1.0}, 'Lipid-body beta-oxidation + glyoxylate'),
    ('04_oleate',        'Oleate (C18:1)',     'pre-infection', {'cpd15269_e0': 1.0}, 'Major lipid-body FA'),
    ('05_hexacosanoate', 'Hexacosanoate (C26)','pre-infection', {'cpd15240_e0': 1.0}, 'Cuticle-wax FA; expected gap in VLCFA b-ox'),
    # ---- Biotrophic / early colonization ----
    ('06_glucose_low',   'Glucose (LOW, -1)',  'biotrophic',    {'cpd00027_e0': 1.0}, 'Apoplastic low glucose; reference'),
    ('07_sucrose_low',   'Sucrose (LOW, -0.5)','biotrophic',    {'cpd00076_e0': 0.5}, 'Phloem sugar; invertase route'),
    ('08_glutamate',     'L-Glutamate',        'biotrophic',    {'cpd00023_e0': 2.0}, 'Apoplastic AA serving C+N'),
    ('09_glutamine',     'L-Glutamine',        'biotrophic',    {'cpd00053_e0': 2.0}, 'Major N-transport AA, 2N + a-KG'),
    ('10_gaba',          'GABA',               'biotrophic',    {'cpd00281_e0': 2.0}, 'Plant defense AA; tests GABA shunt'),
    # ---- Necrotrophic ----
    ('11_galacturonate', 'D-Galacturonate',    'necrotrophic',  {'cpd00280_e0': 5.0}, 'Pectin breakdown; Ashwell pathway'),
    ('12_xylose',        'Xylose',             'necrotrophic',  {'cpd00154_e0': 5.0}, 'Hemicellulose monomer'),
    ('13_arabinose',     'L-Arabinose',        'necrotrophic',  {'cpd00224_e0': 5.0}, 'Arabinoxylan side-chain monomer'),
    ('14_maltose',       'Maltose (starch)',   'necrotrophic',  {'cpd00179_e0': 5.0}, 'Starch breakdown; maltase route'),
    ('15_mixed_aa',      'Mixed AAs + glc',    'necrotrophic',
        {**{aa: 0.5 for aa in AA_EXCHANGES}, 'cpd00027_e0': 0.5}, 'Cytoplasmic release (all 20 AAs + glc)'),
    # ---- Cocktails ----
    ('16_pre_infection_mix', 'Pre-infection mix', 'cocktail',
        {'cpd00794_e0': 0.3, 'cpd00100_e0': 0.3, 'cpd15269_e0': 0.3}, 'Trehalose + glycerol + oleate'),
    ('17_biotrophic_mix', 'Biotrophic mix', 'cocktail',
        {'cpd00027_e0': 0.3, 'cpd00076_e0': 0.2, 'cpd00023_e0': 0.5,
         'cpd00053_e0': 0.5, 'cpd00130_e0': 0.5}, 'glc + suc + Glu + Gln + malate'),
    ('18_necrotrophic_mix', 'Necrotrophic mix', 'cocktail',
        {'cpd00076_e0': 2.0, 'cpd00179_e0': 2.0, 'cpd00280_e0': 1.0,
         'cpd00154_e0': 1.0, 'cpd00224_e0': 0.5,
         **{aa: 0.2 for aa in AA_EXCHANGES}}, 'suc + maltose + galU + xyl + ara + mixed AAs'),
]


# ---------------------------------------------------------------------------

def apply_media(model, c_sources, aerobic):
    """Close every exchange, then open inorganic + condition C-sources + O2."""
    for ex in model.exchanges:
        ex.lower_bound = 0.0
        # NB: do not touch upper bound -- always allow secretion
    for met_id in INORGANIC:
        rid = 'EX_' + met_id
        if rid in [r.id for r in model.exchanges]:
            model.reactions.get_by_id(rid).lower_bound = -1000.0
    # Condition C-sources (uptake rates given as magnitudes -> negate)
    missing = []
    for met_id, rate in c_sources.items():
        rid = 'EX_' + met_id
        if rid in [r.id for r in model.exchanges]:
            model.reactions.get_by_id(rid).lower_bound = -float(rate)
        else:
            missing.append(met_id)
    # O2 condition
    rid = 'EX_' + O2
    if rid in [r.id for r in model.exchanges]:
        model.reactions.get_by_id(rid).lower_bound = -1000.0 if aerobic else 0.0
    return missing


def carbon_atoms_per_cpd(model):
    """C-atom count per metabolite id from formula. Quick parse."""
    out = {}
    for m in model.metabolites:
        f = m.formula or ''
        # match 'C' followed by digits (or end of element), excluding Ca/Cl/Co/Cu/Cd/Cr/Cs
        f2 = re.sub(r'(Ca|Cl|Co|Cu|Cd|Cr|Cs)\d*', '', f)
        match = re.search(r'C(\d*)(?![a-z])', f2)
        if match:
            out[m.id] = int(match.group(1)) if match.group(1) else 1
        else:
            out[m.id] = 0
    return out


def total_C_uptake(model, c_sources, c_atoms):
    """Theoretical max C-uptake (mmol C/gDW/h) from the spec'd uptake rates."""
    total = 0.0
    for met_id, rate in c_sources.items():
        total += rate * c_atoms.get(met_id, 0)
    return total


def actual_C_uptake_from_solution(sol, c_sources, c_atoms):
    """Realized C uptake using flux solution (some sources may not be used)."""
    total = 0.0
    for met_id in c_sources:
        rid = 'EX_' + met_id
        if rid in sol.fluxes.index:
            v = sol.fluxes[rid]
            if v < 0:  # uptake
                total += -v * c_atoms.get(met_id, 0)
    return total


def run_condition(model, cond_id, label, stage, c_sources, notes,
                  c_atoms, aerobic=True):
    """Apply media for one (condition, O2) pair and return a result dict."""
    with model:
        missing = apply_media(model, c_sources, aerobic)
        # bio_gsm as objective
        for r in model.reactions:
            r.objective_coefficient = 1 if r.id == 'bio_gsm' else 0
        sol = model.optimize()
        biomass = sol.objective_value if sol.status == 'optimal' else None
        status  = sol.status

        active = sum(1 for v in sol.fluxes.values if abs(v) > 1e-6) if biomass else 0
        c_in_theor = total_C_uptake(model, c_sources, c_atoms)
        c_in_real  = actual_C_uptake_from_solution(sol, c_sources, c_atoms) if biomass else 0.0
        bio_per_C  = (biomass / c_in_real) if (biomass and c_in_real > 0) else 0.0

        # Per-condition flux dump (only when feasible)
        if biomass and biomass > 1e-6:
            df = pd.DataFrame({'rxn_id': sol.fluxes.index, 'flux': sol.fluxes.values})
            df = df[df['flux'].abs() > 1e-6].sort_values('flux', key=lambda c: c.abs(),
                                                          ascending=False)
            # Attach name / equation / GPR
            def info(rid):
                try:
                    r = model.reactions.get_by_id(rid)
                    return r.name, r.build_reaction_string(use_metabolite_names=True), r.gene_reaction_rule
                except KeyError:
                    return '', '', ''
            df[['name','equation','gpr']] = df['rxn_id'].apply(lambda r: pd.Series(info(r)))
            df.to_csv(f'{PER_DIR}/{cond_id}_{"aerobic" if aerobic else "anaerobic"}.tsv',
                      sep='\t', index=False)

    return {
        'condition_id': cond_id,
        'label': label,
        'stage': stage,
        'O2': 'aerobic' if aerobic else 'anaerobic',
        'c_sources': ';'.join(f'{k}:{v}' for k, v in c_sources.items()),
        'missing_exchanges': ';'.join(missing) if missing else '',
        'status': status,
        'biomass': round(biomass, 6) if biomass is not None else None,
        'active_fluxes': active,
        'C_uptake_theor_mmolC': round(c_in_theor, 3),
        'C_uptake_realized_mmolC': round(c_in_real, 3),
        'biomass_per_C': round(bio_per_C, 6),
        'notes': notes,
    }


def main():
    print(f'loading: {MODEL_PATH}')
    model = cobra.io.load_json_model(MODEL_PATH)
    print(f'  {len(model.reactions)} rxns / {len(model.metabolites)} mets / {len(model.genes)} genes')

    c_atoms = carbon_atoms_per_cpd(model)

    print(f'\nrunning {len(CONDITIONS)} conditions x 2 (aerobic/anaerobic) = {2*len(CONDITIONS)} simulations')
    rows = []
    for cond_id, label, stage, c_sources, notes in CONDITIONS:
        for aerobic in (True, False):
            r = run_condition(model, cond_id, label, stage, c_sources, notes,
                              c_atoms, aerobic=aerobic)
            rows.append(r)
            tag = 'AER' if aerobic else 'ANA'
            bio = r['biomass'] if r['biomass'] is not None else 'n/a'
            print(f'  {cond_id:24s} {tag}: bio={bio!s:>10s}  active={r["active_fluxes"]:>4d}  status={r["status"]}')

    df = pd.DataFrame(rows)
    df.to_csv(SUMMARY, sep='\t', index=False)
    print(f'\nsaved: {SUMMARY}')
    print(f'per-condition flux dumps in: {PER_DIR}/')

    print('\n=== SUMMARY by condition (aerobic biomass) ===')
    print(df[df['O2']=='aerobic'][['condition_id','stage','biomass',
                                     'C_uptake_realized_mmolC','biomass_per_C','status']]
          .to_string(index=False))


if __name__ == '__main__':
    main()
