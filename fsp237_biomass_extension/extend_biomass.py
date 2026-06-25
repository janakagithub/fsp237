#!/usr/bin/env python3
"""Extend the FSP237 ATP-safe GSM biomass with Colletotrichum-specific
melanin, chitin, and ascorbate end-products.

Source spec:  15June2026_Missing reactions GEM_higginsianum.xlsx
              (pathways curated by user from C. higginsianum / CH63R_* genes)

Workflow:
  1. Load fsp237_atp_safe_gsm.json (the curated, ATP-yield-locked model).
  2. Parse the Excel sheet; for each reaction:
       - if id already in model, skip with a note
       - else add the reaction + any missing metabolites, with bounds from the
         'direction' column ('>'  ->  [0, 1000], '<' -> [-1000, 0], '<=>' -> reversible)
  3. Add the one transport gap (D-arabinono-1,4-lactone mito -> cyto) explicitly.
  4. Extend bio_gsm with melanin, chitin, and ascorbate at literature coefficients.
  5. Validate:
       - aerobic + anaerobic biomass with the extended bio_gsm
       - ATP yield (bio1 objective, glc=-1) must stay at 30 / 2
       - new metabolites can actually be produced (per-metabolite demand FBA)
  6. Save extended model + a report tsv.
"""
import csv
import json
import os
import re
import sys

import cobra
import pandas as pd

BASE = '/home/janakae/fungalTemplate/imm904CobraModel'
INPUT_MODEL  = f'{BASE}/fsp237_atp_safe_gsm/fsp237_atp_safe_gsm.json'
GSM_MODEL    = f'{BASE}/fsp237_minimal_glucose/fsp237_minimal_glucose.json'  # for media bounds
XL_PATH      = f'{BASE}/fsp237_biomass_extension/15June2026_Missing reactions GEM_higginsianum.xlsx'
OUT_DIR      = f'{BASE}/fsp237_biomass_extension'
OUTPUT_MODEL = f'{OUT_DIR}/fsp237_atp_safe_gsm_extended.json'
REPORT_TSV   = f'{OUT_DIR}/extension_report.tsv'
BIOMASS_LOG  = f'{OUT_DIR}/biomass_extension_log.tsv'

# ---- biomass coefficients for the new end products --------------------------
# Coefficients are mmol/gDW biomass. Order-of-magnitude estimates from
# fungal cell-wall composition literature; adjust as your wet-lab data warrants.
NEW_BIOMASS = {
    'cpd11683_c0': 0.05,    # Chitin               (~1% w/w cell wall, MW ~203 g/mol)
    'cpd12744_c0': 0.001,   # Melanin              (small but non-zero; melanized cells)
    'cpd00314_c0': 0.20,    # D-Mannitol           (compatible solute; 0.10-0.50 range,
                            #                       starting at 0.20 ~= 3.6% w/w with MW 182)
                            # Solomon, Tan & Oliver 2005 MPMI 18:110-115 (PMID 15720079);
                            # Dijksterhuis et al. 2006 Biochem J (PMID 16987106).
    'cpd12148_c0': 0.05,    # 1,3-alpha-D-Glucan   (major Colletotrichum cell-wall polysaccharide;
                            #                       starts at 0.05 mmol/gDW as initial estimate)
    # cpd03757_c0 (D-erythro-Ascorbate) intentionally NOT in biomass:
    # the C. higginsianum ascorbate pathway in the spec needs D-arabinose
    # (cpd00185) as its starting substrate, which has no producer in the
    # current model and no media uptake. The pathway reactions are added so
    # the route exists if arabinose ever becomes available, but ascorbate is
    # not a biomass demand here.
}

# Rebalance bio_gsm coefficients for metabolites ALREADY in the GSM biomass.
# The loop that consumes NEW_BIOMASS handles override correctly when the
# metabolite already has a coef in bio_gsm; this dict is split out only for
# documentation clarity.
REBALANCED_BIOMASS = {
    # Glycogen: KBase default was 0.5185 mmol/gDW (>10% w/w with MW 162).
    # That overstates intracellular storage; Colletotrichum mycelium typically
    # stores 1-5% w/w glycogen (Solomon 2005 MPMI). 0.30 mmol/gDW ~= 4.9% w/w.
    'cpd00155_c0': 0.30,

    # dNTP rebalance to FSP237 GC content. The KBase default was a textbook
    # 60/40 A+T/G+C split (dAMP=dTMP=0.0036, dGMP=dCMP=0.0024). Per
    # [[project-fsp237-vs-tx430bb]], use TX430BB (Baroncelli 2014, PMID 24926053)
    # as the genome-scale proxy: 52.70% GC -> A+T = 47.30%.
    # With sum = 0.012 mmol/gDW kept constant:
    #   dGMP = dCMP = 0.012 * 0.527 / 2 = 0.003162  -> 0.00316
    #   dAMP = dTMP = 0.012 * 0.473 / 2 = 0.002838  -> 0.00284
    # Sibling-species support: Buiate et al. 2017, BMC Genomics 18:67
    # (PMID 28073340).
    'cpd00294_c0': 0.00284,  # dAMP
    'cpd00298_c0': 0.00284,  # dTMP
    'cpd00296_c0': 0.00316,  # dGMP
    'cpd00206_c0': 0.00316,  # dCMP
}

# Fungal mannitol pathway (gap-fill: not in source GSM or CMM).
# F6P <-> mannitol-1-P <-> mannitol. Two reversible enzymes connect glycolysis
# to mannitol as a compatible solute and NADPH sink.
#
# Fungi (Aspergillus, Colletotrichum, Stagonospora etc.) use the NADP-specific
# mannitol-1-phosphate dehydrogenase (MpdA, EC 1.1.1.138-like), NOT the
# bacterial NAD-version (EC 1.1.1.17, ModelSEED rxn00546 canonical).
# Reference pathway: MetaCyc PWY-3881 (mannitol biosynthesis I).
# Bio role: drains excess cytosolic NADPH from the PPP into mannitol storage.
MANNITOL_RXNS = [
    {
        'rxn_id': 'rxn00546_c0',
        'name': 'D-Mannitol-1-phosphate 5-dehydrogenase (NADP), fungal',
        'direction': '<=>',     # F6P + NADPH + H+ <-> mannitol-1-P + NADP+ (forward in fungi)
        'equation': 'cpd00006_c0 + cpd27436_c0 <=> cpd00005_c0 + cpd00067_c0 + cpd00072_c0',
        #            NADP        mannitol-1-P     NADPH        H+           F6P
        'gpr': '',              # no GPR provided; would be MpdA ortholog in C. sublineola
        'ec': '1.1.1.138',
        'subsystem': 'mannitol biosynthesis (compatible solute, NADPH sink)',
        'kegg': '',
    },
    {
        'rxn_id': 'rxn01560_c0',
        'name': 'D-Mannitol-1-phosphate phosphohydrolase',
        'direction': '<=>',     # MtPP / mannitol-1-phosphatase
        'equation': 'cpd00001_c0 + cpd27436_c0 <=> cpd00009_c0 + cpd00314_c0',
        'gpr': '',
        'ec': '3.1.3.22',
        'subsystem': 'mannitol biosynthesis (compatible solute)',
        'kegg': '',
    },
]

# alpha-1,3-glucan synthase (cell wall polysaccharide; major Colletotrichum component).
# ModelSEED rxn15561 = KEGG R06045, EC 2.4.1.183 (reversible in ModelSEED/KEGG).
#   UDP-glucose (cpd00026) <=> UDP (cpd00014) + 1,3-alpha-D-Glucan (cpd12148)
# Cytosolic (_c0). GPR: FSP237 ortholog of C. higginsianum CH35J_000682
# (Ags1, TID06255.1) -- BLASTp 91.2% id full-length to fsp237 gene_5001
# (qcov 100%, scov 100%, E=0.0, bitscore 4507).
# Biological reference: Fujikawa et al. 2012, PLoS Pathog 8:e1002882
# (PMID 22927818, DOI 10.1371/journal.ppat.1002882).
ALPHA13_GLUCAN_RXNS = [
    {
        'rxn_id': 'rxn15561_c0',
        'name': 'UDP-glucose:alpha-D-(1-3)-glucan 3-alpha-D-glucosyltransferase',
        'direction': '<=>',
        'equation': 'cpd00026_c0 <=> cpd00014_c0 + cpd12148_c0',
        #            UDP-glucose       UDP           1,3-alpha-D-Glucan
        'gpr': 'gene_5001',
        'ec': '2.4.1.183',
        'subsystem': 'cell wall polysaccharide / alpha-1,3-glucan',
        'kegg': 'R06045',
    },
]

# Block the two wrong-substrate "GalNAc" reactions so chitin must flow through
# the user-curated canonical UDP-GlcNAc -> Chitin route (rxn00293 -> rxn15558).
# Both of these reactions are NAMED as if they handle UDP-GlcNAc but actually
# consume/produce cpd00175 (UDP-N-acetyl-D-GALACTOSAMINE), which is the C4
# epimer of UDP-GlcNAc and not a chitin precursor in real biology.
GALNAC_BUG_RXNS = [
    'rxn09486_c0',   # mislabeled "UDP-N-acetylglucosamine pyrophosphorylase";
                     # produces cpd00175 (UDP-GalNAc) instead of cpd00037 (UDP-GlcNAc)
    'rxn09490_c0',   # mislabeled "Chitin synthase 1"; consumes UDP-GalNAc
]

# ---- helpers ----------------------------------------------------------------

DIRECTION_TO_BOUNDS = {
    '>':   (0.0,    1000.0),
    '=>':  (0.0,    1000.0),
    '<':   (-1000.0,    0.0),
    '<=':  (-1000.0,    0.0),
    '<=>': (-1000.0, 1000.0),
    '<>':  (-1000.0, 1000.0),
}

# Equation parser: handles "cpd00007 + (2) cpd00069 => (2) cpd00291"
TERM_RE = re.compile(r'(?:\((\d+(?:\.\d+)?)\)\s*)?(cpd\d+)')
ARROW_RE = re.compile(r'\s*(<=>|=>|<=|->|-->|<-|<--)\s*')


def parse_equation(eq_str, comp_suffix):
    """Return dict {cpd_id_with_comp: signed_coef}."""
    arrow_match = ARROW_RE.search(eq_str)
    if not arrow_match:
        raise ValueError(f'no arrow found in: {eq_str!r}')
    arrow = arrow_match.group(1)
    lhs, rhs = eq_str.split(arrow, 1)
    stoich = {}
    for side, sign in [(lhs, -1), (rhs, +1)]:
        for coef_str, cpd in TERM_RE.findall(side):
            coef = float(coef_str) if coef_str else 1.0
            mid = f'{cpd}{comp_suffix}'
            stoich[mid] = stoich.get(mid, 0.0) + sign * coef
    return stoich


def ensure_metabolite(model, mid, source_model=None):
    """Return the cobra.Metabolite for mid in model, adding it if missing.
    Pulls name/formula/charge/compartment from source_model when possible;
    otherwise infers compartment from id suffix."""
    if mid in [x.id for x in model.metabolites]:
        return model.metabolites.get_by_id(mid)
    comp = mid.rsplit('_', 1)[-1]
    name = mid
    formula = ''
    charge = None
    if source_model is not None and mid in [x.id for x in source_model.metabolites]:
        sm = source_model.metabolites.get_by_id(mid)
        name = sm.name or mid
        formula = sm.formula or ''
        charge = sm.charge
    new_m = cobra.Metabolite(id=mid, name=name, formula=formula,
                              charge=(charge if charge is not None else 0),
                              compartment=comp)
    model.add_metabolites([new_m])
    return new_m


def add_reaction_from_spec(model, rid, name, eq_str, direction, comp_suffix,
                            gpr='', subsystem='', annotation=None,
                            source_model=None):
    """Add a reaction to the model with the given equation.

    If the reaction id already exists in the model, do nothing and return its
    status. Otherwise, build it (creating any missing metabolites first) and
    add it. Returns (status, reason)."""
    if rid in [r.id for r in model.reactions]:
        return ('already_present', '')
    bounds = DIRECTION_TO_BOUNDS.get(direction.strip())
    if bounds is None:
        return ('skipped', f'unknown direction: {direction!r}')
    try:
        stoich = parse_equation(eq_str, comp_suffix)
    except Exception as e:
        return ('skipped', f'parse error: {e}')
    # Create metabolites
    met_dict = {}
    for mid, coef in stoich.items():
        met = ensure_metabolite(model, mid, source_model=source_model)
        met_dict[met] = coef
    rxn = cobra.Reaction(rid, name=name, lower_bound=bounds[0], upper_bound=bounds[1])
    rxn.add_metabolites(met_dict)
    if gpr: rxn.gene_reaction_rule = gpr
    if subsystem: rxn.subsystem = subsystem
    if annotation: rxn.annotation = dict(annotation)
    model.add_reactions([rxn])
    return ('added', '')


def media_atp_yield(model, atp_rxn, glc=1.0, aerobic=True):
    """Standard minimal media for ATP yield test."""
    with model:
        for ex in model.exchanges:
            mid = next(iter(ex.metabolites.keys())).id
            if mid in {'cpd00001_e0','cpd00009_e0','cpd00011_e0','cpd00067_e0',
                       'cpd00013_e0','cpd00048_e0','cpd10515_e0','cpd00971_e0',
                       'cpd00205_e0'}:
                ex.lower_bound = -1000
            else:
                ex.lower_bound = 0
        model.reactions.get_by_id('EX_cpd00027_e0').lower_bound = -glc
        model.reactions.get_by_id('EX_cpd00027_e0').upper_bound = 0
        model.reactions.get_by_id('EX_cpd00007_e0').lower_bound = -1000 if aerobic else 0
        for r in model.reactions: r.objective_coefficient = 1 if r.id == atp_rxn else 0
        return model.optimize().objective_value


def gsm_media_biomass(model, gsm, bio_rxn, aerobic=True):
    """Test biomass on the GSM's own minimal-glucose media bounds."""
    with model:
        gsm_bounds = {r.id: (r.lower_bound, r.upper_bound) for r in gsm.exchanges}
        for ex in model.exchanges:
            if ex.id in gsm_bounds:
                ex.lower_bound, ex.upper_bound = gsm_bounds[ex.id]
            else:
                ex.lower_bound = 0
        if not aerobic:
            model.reactions.get_by_id('EX_cpd00007_e0').lower_bound = 0
        for r in model.reactions: r.objective_coefficient = 1 if r.id == bio_rxn else 0
        return model.optimize().objective_value


def can_produce(model, gsm, met_id, aerobic=True):
    """Max-produce check: add a temporary sink and maximize its flux."""
    with model:
        gsm_bounds = {r.id: (r.lower_bound, r.upper_bound) for r in gsm.exchanges}
        for ex in model.exchanges:
            if ex.id in gsm_bounds:
                ex.lower_bound, ex.upper_bound = gsm_bounds[ex.id]
            else:
                ex.lower_bound = 0
        if not aerobic:
            model.reactions.get_by_id('EX_cpd00007_e0').lower_bound = 0
        if met_id not in [x.id for x in model.metabolites]:
            return 0.0
        met = model.metabolites.get_by_id(met_id)
        sink = cobra.Reaction(f'_DM_{met_id}')
        sink.add_metabolites({met: -1})
        sink.lower_bound = 0; sink.upper_bound = 1000
        model.add_reactions([sink])
        for r in model.reactions: r.objective_coefficient = 1 if r.id == sink.id else 0
        return model.optimize().objective_value or 0.0


# ---- main -------------------------------------------------------------------

def main():
    print(f'loading model: {INPUT_MODEL}')
    model = cobra.io.load_json_model(INPUT_MODEL)
    gsm = cobra.io.load_json_model(GSM_MODEL)
    print(f'  start: {len(model.reactions)} rxns, {len(model.metabolites)} mets')

    # Baseline measurements
    bio_baseline_a = gsm_media_biomass(model, gsm, 'bio_gsm', aerobic=True)
    bio_baseline_n = gsm_media_biomass(model, gsm, 'bio_gsm', aerobic=False)
    atp_baseline_a = media_atp_yield(model, 'bio1', glc=1.0, aerobic=True)
    atp_baseline_n = media_atp_yield(model, 'bio1', glc=1.0, aerobic=False)
    print(f'\nBASELINE biomass (pre-extension): aer={bio_baseline_a:.4f}, ana={bio_baseline_n:.4f}')
    print(f'BASELINE ATP yield                : aer={atp_baseline_a:.3f}, ana={atp_baseline_n:.3f}')

    # ---- Parse Excel ----
    df = pd.read_excel(XL_PATH).dropna(subset=['rxn id'])
    # Drop the placeholder row whose 'rxn id' is literally the string "rxn"
    df = df[df['rxn id'].astype(str).str.strip() != 'rxn']

    report_rows = []
    print(f'\n=== adding {len(df)} reactions from spec ===')
    for _, row in df.iterrows():
        rid = str(row['rxn id']).strip()
        comp_suffix = '_' + str(row['compartment']).strip().replace('_', '')
        # Excel has 'c_0'/'m_0'; normalize to '_c0'/'_m0'
        comp_suffix = '_' + comp_suffix.replace('__', '_').lstrip('_').replace('_', '', 0)
        # Simpler: just take the suffix from rid itself when available
        suffix_match = re.search(r'_([a-z]\d+)$', rid)
        if suffix_match:
            comp_suffix = '_' + suffix_match.group(1)
        else:
            comp_suffix = '_c0'

        name = str(row['name (ModelSEED)']).strip() if pd.notna(row['name (ModelSEED)']) else rid
        eq   = str(row['equation']).strip() if pd.notna(row['equation']) else ''
        direc = str(row['direction']).strip() if pd.notna(row['direction']) else '>'
        gpr_raw = row['gene ID KEGG;  Uniprot; NCBI']
        gpr = ''
        if pd.notna(gpr_raw):
            # collapse whitespace; semicolons -> ' or '
            tokens = [t.strip() for t in re.split(r'[;\n\t]+', str(gpr_raw)) if t.strip()]
            gpr = ' or '.join(tokens)
        ec = str(row['enzyme EC']).strip() if pd.notna(row['enzyme EC']) else ''
        subsystem = str(row['subsystem']).strip() if pd.notna(row['subsystem']) else ''
        kegg = str(row['kegg id']).strip() if pd.notna(row['kegg id']) else ''
        annotation = {}
        if kegg: annotation['kegg.reaction'] = kegg
        if ec: annotation['ec-code'] = ec

        status, reason = add_reaction_from_spec(
            model, rid=rid, name=name, eq_str=eq, direction=direc,
            comp_suffix=comp_suffix, gpr=gpr, subsystem=subsystem,
            annotation=annotation, source_model=gsm,
        )
        report_rows.append({
            'rxn_id': rid, 'status': status, 'reason': reason,
            'direction': direc, 'subsystem': subsystem,
            'name': name, 'equation': eq, 'gpr': gpr, 'ec': ec,
        })
        flag = {'added': '+', 'already_present': '·', 'skipped': '!'}[status]
        print(f'  {flag} {rid:18s} [{status:15s}] {name[:60]}')
        if reason: print(f'      reason: {reason}')

    # ---- Add the fungal mannitol pathway (hardcoded; not in Excel) ----
    print(f'\n=== adding fungal mannitol pathway ({len(MANNITOL_RXNS)} rxns) ===')
    for spec in MANNITOL_RXNS:
        rid = spec['rxn_id']
        annotation = {}
        if spec.get('kegg'): annotation['kegg.reaction'] = spec['kegg']
        if spec.get('ec'):   annotation['ec-code']      = spec['ec']
        status, reason = add_reaction_from_spec(
            model, rid=rid, name=spec['name'], eq_str=spec['equation'],
            direction=spec['direction'], comp_suffix='_c0',
            gpr=spec.get('gpr', ''), subsystem=spec.get('subsystem', ''),
            annotation=annotation, source_model=gsm,
        )
        report_rows.append({
            'rxn_id': rid, 'status': status, 'reason': reason or 'mannitol pathway (Solomon 2005, Dijksterhuis 2006)',
            'direction': spec['direction'], 'subsystem': spec.get('subsystem', ''),
            'name': spec['name'], 'equation': spec['equation'],
            'gpr': spec.get('gpr', ''), 'ec': spec.get('ec', ''),
        })
        flag = {'added': '+', 'already_present': '·', 'skipped': '!'}[status]
        print(f'  {flag} {rid:18s} [{status:15s}] {spec["name"][:60]}')
        if reason: print(f'      reason: {reason}')

    # ---- Add the alpha-1,3-glucan synthase (hardcoded; CH35J_000682 ortholog) ----
    print(f'\n=== adding alpha-1,3-glucan synthase ({len(ALPHA13_GLUCAN_RXNS)} rxn) ===')
    for spec in ALPHA13_GLUCAN_RXNS:
        rid = spec['rxn_id']
        annotation = {}
        if spec.get('kegg'): annotation['kegg.reaction'] = spec['kegg']
        if spec.get('ec'):   annotation['ec-code']      = spec['ec']
        status, reason = add_reaction_from_spec(
            model, rid=rid, name=spec['name'], eq_str=spec['equation'],
            direction=spec['direction'], comp_suffix='_c0',
            gpr=spec.get('gpr', ''), subsystem=spec.get('subsystem', ''),
            annotation=annotation, source_model=gsm,
        )
        # cpd12148 (1,3-alpha-D-Glucan) is not in the source GSM, so ensure_metabolite
        # creates it as a bare id with no name/formula. Patch from ModelSEED (KEGG C02616).
        if 'cpd12148_c0' in [x.id for x in model.metabolites]:
            polymer = model.metabolites.get_by_id('cpd12148_c0')
            if not polymer.name or polymer.name == polymer.id:
                polymer.name = '1,3-alpha-D-Glucan_c0'
            if not polymer.formula:
                polymer.formula = 'C6H10O5R2'   # per glucose unit; polymer placeholder R groups
            polymer.annotation = {'kegg.compound': 'C02616', 'modelseed.compound': 'cpd12148'}
        report_rows.append({
            'rxn_id': rid, 'status': status,
            'reason': reason or 'alpha-1,3-glucan synthase; GPR=gene_5001 (91.2% id to CH35J_000682 Ags1); Fujikawa 2012 PMID 22927818',
            'direction': spec['direction'], 'subsystem': spec.get('subsystem', ''),
            'name': spec['name'], 'equation': spec['equation'],
            'gpr': spec.get('gpr', ''), 'ec': spec.get('ec', ''),
        })
        flag = {'added': '+', 'already_present': '·', 'skipped': '!'}[status]
        print(f'  {flag} {rid:18s} [{status:15s}] {spec["name"][:60]}')
        if reason: print(f'      reason: {reason}')

    # ---- Add the missing transport: D-arabinono-1,4-lactone mito -> cyto ----
    # cpd00496 = D-Arabinono-1,4-lactone
    transport_id = 'tx_cpd00496_mc'
    if transport_id not in [r.id for r in model.reactions]:
        m_m = ensure_metabolite(model, 'cpd00496_m0', source_model=gsm)
        m_c = ensure_metabolite(model, 'cpd00496_c0', source_model=gsm)
        r = cobra.Reaction(transport_id,
                            name='D-Arabinono-1,4-lactone transport (mito -> cyto)',
                            lower_bound=-1000.0, upper_bound=1000.0)
        r.add_metabolites({m_m: -1, m_c: 1})
        model.add_reactions([r])
        report_rows.append({
            'rxn_id': transport_id, 'status': 'added', 'reason': 'transport stub from spec',
            'direction': '<=>', 'subsystem': 'Ascorbate transport',
            'name': r.name, 'equation': 'cpd00496_m0 <=> cpd00496_c0',
            'gpr': '', 'ec': '',
        })
        print(f'  + {transport_id} [added (manual transport stub)]')

    print(f'\n  model now: {len(model.reactions)} rxns, {len(model.metabolites)} mets')

    # ---- Block the GalNAc-bug reactions so chitin flows through canonical pathway ----
    print(f'\n=== blocking GalNAc-bug reactions (forces canonical chitin pathway) ===')
    for rid in GALNAC_BUG_RXNS:
        if rid in [r.id for r in model.reactions]:
            r = model.reactions.get_by_id(rid)
            old = r.bounds
            r.bounds = (0.0, 0.0)
            print(f'  {rid}: bounds {old} -> (0.0, 0.0)  [{r.name}]')
            report_rows.append({
                'rxn_id': rid, 'status': 'blocked', 'reason': 'wrong-substrate (UDP-GalNAc instead of UDP-GlcNAc)',
                'direction': 'blocked', 'subsystem': 'chitin synthesis (corrected)',
                'name': r.name, 'equation': r.build_reaction_string(use_metabolite_names=False),
                'gpr': r.gene_reaction_rule, 'ec': '',
            })

    # ---- Extend / rebalance bio_gsm coefficients ----
    print(f'\n=== extending bio_gsm with new biomass components ===')
    bio = model.reactions.get_by_id('bio_gsm')
    bio_log = []

    def _apply_biomass(mid, coef, kind):
        if mid not in [x.id for x in model.metabolites]:
            print(f'  ! {mid} not in model after additions; SKIPPING (check pathway)')
            bio_log.append({'metabolite': mid, 'coefficient': 0, 'status': f'skipped (missing, {kind})'})
            return
        met = model.metabolites.get_by_id(mid)
        existing = bio.metabolites.get(met, 0)
        if existing != 0:
            print(f'  · {mid} already in bio_gsm at coef {existing:.5f}; updating to {-coef:.5f} ({kind})')
            bio.add_metabolites({met: -coef - existing})  # net = -coef
        else:
            bio.add_metabolites({met: -coef})
            print(f'  + {mid} added to bio_gsm: -{coef} ({met.name}) ({kind})')
        bio_log.append({'metabolite': mid, 'name': met.name,
                         'coefficient': -coef, 'status': kind})

    for mid, coef in NEW_BIOMASS.items():
        _apply_biomass(mid, coef, kind='added')

    print(f'\n=== rebalancing existing bio_gsm coefficients ===')
    for mid, coef in REBALANCED_BIOMASS.items():
        _apply_biomass(mid, coef, kind='rebalanced')

    # ---- Validate ----
    print('\n=== POST-EXTENSION VALIDATION ===')
    bio_a = gsm_media_biomass(model, gsm, 'bio_gsm', aerobic=True)
    bio_n = gsm_media_biomass(model, gsm, 'bio_gsm', aerobic=False)
    atp_a = media_atp_yield(model, 'bio1', glc=1.0, aerobic=True)
    atp_n = media_atp_yield(model, 'bio1', glc=1.0, aerobic=False)
    print(f'biomass  aer={bio_a:.4f}  (was {bio_baseline_a:.4f})')
    print(f'biomass  ana={bio_n:.4f}  (was {bio_baseline_n:.4f})')
    print(f'ATP/glc  aer={atp_a:.3f}  (target 30)')
    print(f'ATP/glc  ana={atp_n:.3f}  (target 2)')

    print('\nCan each new/rebalanced biomass component be produced?')
    for mid in list(NEW_BIOMASS) + list(REBALANCED_BIOMASS):
        prod_a = can_produce(model, gsm, mid, aerobic=True)
        prod_n = can_produce(model, gsm, mid, aerobic=False)
        nm = model.metabolites.get_by_id(mid).name if mid in [x.id for x in model.metabolites] else '(missing)'
        print(f'  {mid:15s} {nm[:30]:30s}  max-prod aer={prod_a:.4f}, ana={prod_n:.4f}')

    # Save
    cobra.io.save_json_model(model, OUTPUT_MODEL)
    print(f'\nsaved: {OUTPUT_MODEL}')
    pd.DataFrame(report_rows).to_csv(REPORT_TSV, sep='\t', index=False)
    print(f'saved: {REPORT_TSV}')
    pd.DataFrame(bio_log).to_csv(BIOMASS_LOG, sep='\t', index=False)
    print(f'saved: {BIOMASS_LOG}')

    # Final summary block
    print('\n=== SUMMARY ===')
    added = sum(1 for r in report_rows if r['status'] == 'added')
    present = sum(1 for r in report_rows if r['status'] == 'already_present')
    skipped = sum(1 for r in report_rows if r['status'] == 'skipped')
    print(f'  reactions added           : {added}')
    print(f'  reactions already present : {present}')
    print(f'  reactions skipped         : {skipped}')
    print(f'  biomass aerobic preserved : {abs(bio_a - bio_baseline_a) < 0.01} (delta {bio_a-bio_baseline_a:+.4f})')
    print(f'  ATP yield 30/2 preserved  : aer={abs(atp_a-30)<0.01} ana={abs(atp_n-2)<0.01}')


if __name__ == '__main__':
    main()
