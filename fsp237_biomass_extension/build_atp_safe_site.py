#!/usr/bin/env python3
"""Regenerate /home/janakae/fsp237/atp-safe/reactions.json from the current
extended GSM. Preserves the existing schema (consumed by atp-safe/index.html
via DataTables).

Schema produced:
  top-level: generated_at, model_id, model_name, n_reactions,
             n_active_aerobic, n_active_anaerobic,
             biomass_aerobic, biomass_anaerobic,
             atp_yield_aerobic, atp_yield_anaerobic,
             reactions: [...], biomass: {...}
  per-reaction: id, name, compartment, compartment_name, direction,
                lb, ub, equation_cpd, equation_name, ec, gpr, subsystem,
                flux_aerobic, flux_anaerobic
  biomass: { reaction_id, components: [{id, name, formula, compartment,
            coefficient, role, category}] }
"""
from datetime import datetime, timezone
import json
import os
import sys

import cobra

BASE = '/home/janakae/fungalTemplate/imm904CobraModel'
# Source the LATEST curated model: V6 = gapfilled + gene-integrated + deduped
# + direction-locked (degradation-only β-ox + Ashwell). The GPR overhaul
# from gpr-update/ is already inside this -- V6 supersedes the
# fsp237_atp_safe_gsm_gpr_updated.json file. See
# simulations/gapfill_v1_v2/reports/SUMMARY.md for the version trail.
MODEL_PATH = f'{BASE}/simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json'
GSM_PATH   = f'{BASE}/fsp237_minimal_glucose/fsp237_minimal_glucose.json'
OUT_PATH   = '/home/janakae/fsp237/atp-safe/reactions.json'

# Simulation-panel results to embed under reactions.json -> 'simulations'.
# Optional -- the build still works if these files are missing (the
# Simulations tab on the site will then show 'no data').
SIM_TSV    = f'{BASE}/simulations/simulation_results.tsv'
# Per-condition literature support (compiled in condition_literature.tsv).
# Columns: condition_id, pmid, short_citation, key_finding.
# A condition may have multiple rows (multiple citations).
LIT_TSV    = f'{BASE}/simulations/condition_literature.tsv'

# Keep these stable so the static site title/badge stays consistent.
MODEL_ID   = 'Colletotrichum_sublineola_FSP237_MyCoCosm_draftModel'
MODEL_NAME = 'Colletotrichum sublineola'

COMPARTMENT_NAMES = {
    'c0': 'cytosol', 'c': 'cytosol',
    'm0': 'mitochondria', 'm': 'mitochondria',
    'e0': 'extracellular', 'e': 'extracellular',
    'r0': 'endoplasmic reticulum',
    'g0': 'Golgi',
    'n0': 'nucleus',
    'x0': 'peroxisome',
    'v0': 'vacuole',
}


# Curated category map for the well-known biomass components. Anything not
# listed here falls through to keyword/formula heuristics in categorize().
BIOMASS_CATEGORY_BY_ID = {
    # Energy / growth-associated maintenance (GAM stoich block)
    'cpd00001': 'Energy / GAM',  # H2O
    'cpd00002': 'Energy / GAM',  # ATP
    'cpd00008': 'Energy / GAM',  # ADP
    'cpd00009': 'Energy / GAM',  # Phosphate
    'cpd00067': 'Energy / GAM',  # H+
    # RNA mononucleotides
    'cpd00018': 'RNA', 'cpd00046': 'RNA', 'cpd00091': 'RNA', 'cpd00126': 'RNA',
    # DNA mononucleotides
    'cpd00206': 'DNA', 'cpd00294': 'DNA', 'cpd00296': 'DNA', 'cpd00298': 'DNA',
    # Cell wall polysaccharides
    'cpd11683': 'Cell wall',   # Chitin
    'cpd11685': 'Cell wall',   # Mannan
    'cpd11791': 'Cell wall',   # 1,3-beta-Glucan
    'cpd12148': 'Cell wall',   # 1,3-alpha-D-Glucan
    # Storage carbohydrates
    'cpd00155': 'Storage',     # Glycogen
    'cpd00794': 'Storage',     # Trehalose
    'cpd00314': 'Storage',     # D-Mannitol (compatible solute)
    # Pigments / other Colletotrichum-specific
    'cpd12744': 'Pigment',     # Melanin
    # Sterols (membrane)
    'cpd01170': 'Sterol',      # Ergosterol
    'cpd03221': 'Sterol',      # Zymosterol
    'cpd14514': 'Sterol',      # Episterol
    # Cofactors / vitamins
    'cpd00003': 'Cofactor', 'cpd00006': 'Cofactor', 'cpd00010': 'Cofactor',
    'cpd00015': 'Cofactor', 'cpd00016': 'Cofactor', 'cpd00017': 'Cofactor',
    'cpd00056': 'Cofactor', 'cpd00087': 'Cofactor', 'cpd00104': 'Cofactor',
    'cpd00220': 'Cofactor', 'cpd00345': 'Cofactor', 'cpd00557': 'Cofactor',
    'cpd00644': 'Cofactor', 'cpd11312': 'Cofactor', 'cpd15290': 'Cofactor',
    # Inorganic
    'cpd00048': 'Inorganic',   # Sulfate
    # Membrane phospholipid biomass placeholders (yeast-derived; "_BS" tag)
    'cpd29687': 'Lipid / membrane',  # ps_BS  (phosphatidylserine, biomass-specific)
    'cpd29688': 'Lipid / membrane',  # psetha_BS (phosphatidylethanolamine, biomass-specific)
}

AMINO_ACID_IDS = {
    'cpd00023','cpd00033','cpd00035','cpd00039','cpd00041','cpd00051',
    'cpd00053','cpd00054','cpd00060','cpd00065','cpd00066','cpd00069',
    'cpd00084','cpd00107','cpd00119','cpd00129','cpd00132','cpd00156',
    'cpd00161','cpd00322',
}

CATEGORY_ORDER = [
    'Energy / GAM', 'Amino acid', 'RNA', 'DNA', 'Cell wall', 'Storage',
    'Lipid / membrane', 'Sterol', 'Cofactor', 'Pigment', 'Inorganic', 'Other',
]


def categorize(met):
    """Best-effort category for a biomass metabolite."""
    base = met.id.split('_')[0]
    if base in BIOMASS_CATEGORY_BY_ID:
        return BIOMASS_CATEGORY_BY_ID[base]
    if base in AMINO_ACID_IDS:
        return 'Amino acid'
    nm = (met.name or '').lower()
    formula = (met.formula or '')
    if any(k in nm for k in ('ceramide', 'lecithin', 'phosphatidyl', 'triglyceride',
                              'lipid', 'cardiolipin', 'phosphatidate',
                              'phospholipid', 'phosphatidylserine')):
        return 'Lipid / membrane'
    if any(k in nm for k in ('sterol', 'sterone', 'ergost')):
        return 'Sterol'
    if 'glucan' in nm or 'chitin' in nm or 'mannan' in nm:
        return 'Cell wall'
    if any(k in nm for k in ('coa', 'nad', 'fad', 'biotin', 'thiamin', 'folate',
                              'flavin', 'sam', 'pyridoxal', 'heme', 'ubiquinone',
                              'siroheme', 'phosphopantetheine', 'tetrahydrof')):
        return 'Cofactor'
    if any(k in nm for k in ('mp_', 'monophosphate')) and 'amp' in nm:
        return 'RNA'
    if formula in {'HO4P', 'H2O', 'H', 'O4S'}:
        return 'Inorganic'
    return 'Other'


def direction_from_bounds(lb, ub):
    if lb < 0 and ub > 0: return '<=>'
    if ub > 0 and lb >= 0: return '=>'
    if lb < 0 and ub <= 0: return '<='
    return '='


def primary_compartment(rxn):
    comps = sorted({m.id.rsplit('_', 1)[-1] for m in rxn.metabolites if '_' in m.id})
    return comps[0] if comps else ''


def fba_with_gsm_media(model, gsm, objective_id, aerobic=True):
    """Set GSM media and aerobic/anaerobic O2, optimize `objective_id`,
    return (objective_value, fluxes_dict)."""
    with model:
        gsm_bounds = {r.id: (r.lower_bound, r.upper_bound) for r in gsm.exchanges}
        for ex in model.exchanges:
            if ex.id in gsm_bounds:
                ex.lower_bound, ex.upper_bound = gsm_bounds[ex.id]
            else:
                ex.lower_bound = 0
        if not aerobic and 'EX_cpd00007_e0' in [r.id for r in model.reactions]:
            model.reactions.get_by_id('EX_cpd00007_e0').lower_bound = 0
        for r in model.reactions:
            r.objective_coefficient = 1 if r.id == objective_id else 0
        sol = model.optimize()
        return sol.objective_value, dict(sol.fluxes)


def atp_yield(model, atp_rxn='bio1', aerobic=True, glc=1.0):
    """Minimal media ATP yield per mmol glucose (matches extend_biomass logic)."""
    with model:
        for ex in model.exchanges:
            mid = next(iter(ex.metabolites.keys())).id
            if mid in {'cpd00001_e0','cpd00009_e0','cpd00011_e0','cpd00067_e0',
                       'cpd00013_e0','cpd00048_e0','cpd10515_e0','cpd00971_e0',
                       'cpd00205_e0'}:
                ex.lower_bound = -1000
            else:
                ex.lower_bound = 0
        if 'EX_cpd00027_e0' in [r.id for r in model.reactions]:
            model.reactions.get_by_id('EX_cpd00027_e0').lower_bound = -glc
            model.reactions.get_by_id('EX_cpd00027_e0').upper_bound = 0
        if 'EX_cpd00007_e0' in [r.id for r in model.reactions]:
            model.reactions.get_by_id('EX_cpd00007_e0').lower_bound = -1000 if aerobic else 0
        for r in model.reactions:
            r.objective_coefficient = 1 if r.id == atp_rxn else 0
        return model.optimize().objective_value or 0.0


def main():
    print(f'loading model: {MODEL_PATH}')
    model = cobra.io.load_json_model(MODEL_PATH)
    print(f'loading GSM (for media): {GSM_PATH}')
    gsm = cobra.io.load_json_model(GSM_PATH)

    print('running aerobic FBA on bio_gsm ...')
    bio_a, flx_a = fba_with_gsm_media(model, gsm, 'bio_gsm', aerobic=True)
    print(f'  biomass aerobic   = {bio_a:.4f}')
    print('running anaerobic FBA on bio_gsm ...')
    bio_n, flx_n = fba_with_gsm_media(model, gsm, 'bio_gsm', aerobic=False)
    print(f'  biomass anaerobic = {bio_n:.4f}')

    print('computing ATP yields (bio1, glc=-1) ...')
    atp_a = atp_yield(model, 'bio1', aerobic=True,  glc=1.0)
    atp_n = atp_yield(model, 'bio1', aerobic=False, glc=1.0)
    print(f'  ATP/glc aerobic   = {atp_a:.3f}')
    print(f'  ATP/glc anaerobic = {atp_n:.3f}')

    n_active_a = sum(1 for v in flx_a.values() if abs(v) > 1e-6)
    n_active_n = sum(1 for v in flx_n.values() if abs(v) > 1e-6)

    rxn_rows = []
    for r in model.reactions:
        comp = primary_compartment(r)
        rxn_rows.append({
            'id': r.id,
            'name': r.name,
            'compartment': comp,
            'compartment_name': COMPARTMENT_NAMES.get(comp, comp),
            'direction': direction_from_bounds(r.lower_bound, r.upper_bound),
            'lb': float(r.lower_bound),
            'ub': float(r.upper_bound),
            'equation_cpd': r.build_reaction_string(use_metabolite_names=False),
            'equation_name': r.build_reaction_string(use_metabolite_names=True),
            'ec': r.annotation.get('ec-code', '') if isinstance(r.annotation, dict) else '',
            'gpr': r.gene_reaction_rule,
            'subsystem': r.subsystem or '',
            'flux_aerobic':   float(flx_a.get(r.id, 0.0)),
            'flux_anaerobic': float(flx_n.get(r.id, 0.0)),
        })

    # Extract biomass composition (bio_gsm is the real biomass; bio1 is the
    # ATP-demand objective and is intentionally excluded).
    #
    # Some metabolites came in with no display name or with a redundant
    # "_c0" suffix from the biomass-extension pipeline (ensure_metabolite
    # fallback). Patch them here so the Biomass tab renders human-readable
    # names. The model JSON is untouched -- this is a presentation-layer
    # fix only.
    NAME_OVERRIDES = {
        'cpd00314_c0': 'D-Mannitol',
        'cpd12744_c0': 'Melanin',
        'cpd29687_c0': 'Phosphatidylserine (biomass-spec.)',
        'cpd29688_c0': 'Phosphatidylethanolamine (biomass-spec.)',
        'cpd00794_c0': 'Trehalose',
    }
    # Mark biomass-extension compounds with a visual flag so users can spot
    # them quickly (per fsp237_biomass_extension/biomass_extension_log.tsv).
    EXTENSION_COMPOUNDS = {
        'cpd11683_c0', 'cpd12744_c0', 'cpd00314_c0', 'cpd12148_c0',
        'cpd00155_c0', 'cpd00206_c0', 'cpd00294_c0', 'cpd00296_c0',
        'cpd00298_c0',
    }

    def clean_name(met):
        name = NAME_OVERRIDES.get(met.id) or (met.name or met.id)
        # strip a redundant "_c0" / "_m0" suffix when the compartment is
        # already shown in its own column
        for suf in (f'_{met.compartment}',):
            if suf and name.endswith(suf):
                name = name[: -len(suf)]
                break
        return name

    bio_rxn = model.reactions.get_by_id('bio_gsm')
    biomass_components = []
    for met, coef in sorted(bio_rxn.metabolites.items(), key=lambda kv: kv[0].id):
        biomass_components.append({
            'id': met.id,
            'name': clean_name(met),
            'formula': met.formula or '',
            'compartment': met.compartment or '',
            'coefficient': float(coef),
            'role': 'substrate' if coef < 0 else 'product',
            'category': categorize(met),
            'extension': met.id in EXTENSION_COMPOUNDS,
        })

    # Optional simulation-panel section (from simulations/simulation_results.tsv)
    sim_payload = None
    if os.path.exists(SIM_TSV):
        # Build cpd_id -> human-readable name lookup from the model, so the
        # site can show "cpd00076_e0:2.0 (Sucrose)" instead of just the id.
        name_lookup = {}
        for met in model.metabolites:
            base_name = (met.name or met.id)
            # strip the trailing "_e0" / "_c0" tag if present on the name
            for suf in (f'_{met.compartment}',):
                if base_name.endswith(suf):
                    base_name = base_name[: -len(suf)]
            name_lookup[met.id] = base_name

        def expand(c_sources_str):
            """'cpd00076_e0:2.0;cpd00179_e0:2.0' ->
                'cpd00076_e0:2.0 (Sucrose);cpd00179_e0:2.0 (Maltose)' """
            if not c_sources_str:
                return ''
            out = []
            for term in c_sources_str.split(';'):
                term = term.strip()
                if not term:
                    continue
                if ':' in term:
                    mid, rate = term.split(':', 1)
                    name = name_lookup.get(mid, '')
                    out.append(f'{mid}:{rate} ({name})' if name else term)
                else:
                    out.append(term)
            return ';'.join(out)

        # Read 36 (condition x O2) rows from TSV, then PIVOT to one row per
        # condition with aerobic + anaerobic side-by-side. The site renders
        # condition-rows so aer/ana growth is easy to compare at a glance.
        raw_rows = []
        with open(SIM_TSV) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            for line in fh:
                if not line.strip():
                    continue
                vals = line.rstrip('\n').split('\t')
                row = dict(zip(header, vals))
                for k in ('biomass','C_uptake_theor_mmolC','C_uptake_realized_mmolC',
                          'biomass_per_C'):
                    if k in row and row[k] not in ('', 'None'):
                        try: row[k] = float(row[k])
                        except ValueError: pass
                if 'active_fluxes' in row and row['active_fluxes']:
                    try: row['active_fluxes'] = int(row['active_fluxes'])
                    except ValueError: pass
                if 'c_sources' in row:
                    row['c_sources_pretty'] = expand(row['c_sources'])
                raw_rows.append(row)

        # Load per-condition literature support (1+ citations per condition)
        lit_by_cond = {}
        if os.path.exists(LIT_TSV):
            with open(LIT_TSV) as lfh:
                lh = lfh.readline().rstrip('\n').split('\t')
                for line in lfh:
                    parts = line.rstrip('\n').split('\t')
                    if not parts or not parts[0]: continue
                    lrow = dict(zip(lh, parts))
                    cid = lrow.get('condition_id')
                    if not cid: continue
                    lit_by_cond.setdefault(cid, []).append({
                        'pmid':           lrow.get('pmid', ''),
                        'short_citation': lrow.get('short_citation', ''),
                        'key_finding':    lrow.get('key_finding', ''),
                    })
            print(f'  loaded {sum(len(v) for v in lit_by_cond.values())} citations '
                   f'across {len(lit_by_cond)} conditions from {LIT_TSV}')

        # Group by condition_id, attach _aer and _ana variants.
        by_cond = {}
        for r in raw_rows:
            cid = r['condition_id']
            slot = 'aer' if r.get('O2') == 'aerobic' else 'ana'
            if cid not in by_cond:
                by_cond[cid] = {
                    'condition_id': cid,
                    'label': r.get('label',''),
                    'stage': r.get('stage',''),
                    'c_sources': r.get('c_sources',''),
                    'c_sources_pretty': r.get('c_sources_pretty',''),
                    'notes': r.get('notes',''),
                    'literature': lit_by_cond.get(cid, []),
                }
            entry = by_cond[cid]
            entry[f'biomass_{slot}']        = r.get('biomass')
            entry[f'active_{slot}']         = r.get('active_fluxes')
            entry[f'status_{slot}']         = r.get('status','')
            entry[f'C_uptake_{slot}']       = r.get('C_uptake_realized_mmolC')
            entry[f'biomass_per_C_{slot}']  = r.get('biomass_per_C')
            entry[f'C_uptake_theor_{slot}'] = r.get('C_uptake_theor_mmolC')

        # Preserve condition_id order from the original TSV
        pivot_rows = []
        seen = set()
        for r in raw_rows:
            cid = r['condition_id']
            if cid not in seen:
                pivot_rows.append(by_cond[cid])
                seen.add(cid)

        sim_payload = {
            'source': 'simulations/run_simulation_panel.py',
            'plan': 'fsp237_biomass_extension/INFECTION_SIM_PLAN.md',
            'n_conditions': len(pivot_rows),
            'n_simulations': len(raw_rows),
            'rows': pivot_rows,  # pivoted: one row per condition, aer + ana side-by-side
        }
        print(f'embedded simulations: {len(pivot_rows)} conditions ({len(raw_rows)} simulations) from {SIM_TSV}')

    payload = {
        'generated_at': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        'model_id': MODEL_ID,
        'model_name': MODEL_NAME,
        'n_reactions': len(model.reactions),
        'n_active_aerobic': n_active_a,
        'n_active_anaerobic': n_active_n,
        'biomass_aerobic': float(bio_a),
        'biomass_anaerobic': float(bio_n),
        'atp_yield_aerobic': float(atp_a),
        'atp_yield_anaerobic': float(atp_n),
        'reactions': rxn_rows,
        'biomass': {
            'reaction_id': 'bio_gsm',
            'lb': float(bio_rxn.lower_bound),
            'ub': float(bio_rxn.upper_bound),
            'n_components': len(biomass_components),
            'category_order': CATEGORY_ORDER,
            'components': biomass_components,
        },
        'simulations': sim_payload,
    }

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    # minified to match the format previously committed to the site repo
    with open(OUT_PATH, 'w') as fh:
        json.dump(payload, fh, separators=(',', ':'))
    print(f'\nwrote: {OUT_PATH}')
    print(f'  reactions     : {payload["n_reactions"]}')
    print(f'  active aerobic: {payload["n_active_aerobic"]}')
    print(f'  active anaer. : {payload["n_active_anaerobic"]}')


if __name__ == '__main__':
    main()
