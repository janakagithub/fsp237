#!/usr/bin/env python3
"""Upload V6 (gapfilled + dedup + dirlock + gene-integrated) FSP237 model
to KBase workspace 28277 ("Fungal Biomass testing Narrative").

We do a custom convert that fixes two issues with the stock
cobrakbase.convert_to_kbase:
  1. cobrakbase only declares 'c' / 'e' compartments; FSP237 uses
     c0 / e0 / m0 / x0 / r0 / g0 / n0 / v0 -- we declare all of them.
  2. cobrakbase emits biomasses=[] -- we detect biomass-like reactions
     (`bio_gsm`, `bio1`) and emit them as proper KBase biomass entries.

Requires ~/.kbase/token (auto-loaded by cobrakbase.KBaseAPI()).
"""
import os
import re

import cobra
import cobrakbase
from cobrakbase.core.cobra_to_kbase import (
    convert_to_kbase_reaction, get_compound_references,
)

BASE     = '/home/janakae/fungalTemplate/imm904CobraModel'
MODEL_PATH = f'{BASE}/simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json'
# (Older snapshots are preserved on disk; see SUMMARY.md for the version
#  trail. This script uploads whichever path is set above as MODEL_PATH.)

WS_ID       = 28277          # janakakbase:narrative_1518190880851
OBJECT_ID   = 'fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated'
OBJECT_TYPE = 'KBaseFBA.FBAModel'

# Compartments present in the FSP237 model
COMPARTMENT_SPEC = {
    'c0': ('cytosol',           'Cytosol'),
    'e0': ('extracellular',     'Extracellular'),
    'm0': ('mitochondria',      'Mitochondria'),
    'x0': ('peroxisome',        'Peroxisome'),
    'r0': ('endoplasmic reticulum', 'Endoplasmic reticulum'),
    'g0': ('Golgi',             'Golgi'),
    'n0': ('nucleus',           'Nucleus'),
    'v0': ('vacuole',           'Vacuole'),
}

# Reactions to surface as KBase biomass objects (degradation excluded).
BIOMASS_RXNS = {
    'bio_gsm': 'biomass [from GSM]',
    'bio1':    'ATP demand (legacy bio1 objective)',
}


def build_kbase_model(object_id, m):
    compounds_to_refs = get_compound_references(m)

    # --- modelcompartments ----------------------------------------------------
    used = {met.compartment for met in m.metabolites if met.compartment}
    modelcompartments = []
    for cid in sorted(used):
        name, label = COMPARTMENT_SPEC.get(cid, (cid, cid))
        modelcompartments.append({
            'compartment_ref': f'~/template/compartments/id/{cid[0]}',
            'compartmentIndex': int(re.findall(r'\d+', cid)[0]) if re.search(r'\d+', cid) else 0,
            'id': cid,
            'label': label,
            'pH': 7,
            'potential': 0,
        })

    compartment_to_ref = {cid: f'~/modelcompartments/id/{cid}' for cid in used}

    # --- modelcompounds (with proper per-met compartment refs) ---------------
    modelcompounds = []
    for met in m.metabolites:
        cpd_base = met.id.split('_')[0] if met.id.startswith('cpd') else 'cpd00000'
        comp_ref = compartment_to_ref.get(met.compartment, '~/modelcompartments/id/c0')
        modelcompounds.append({
            'aliases': [],
            'charge': met.charge if met.charge is not None else 0,
            'compound_ref': f'~/template/compounds/id/{cpd_base}',
            'dblinks': {},
            'formula': met.formula or '*',
            'id': met.id,
            'modelcompartment_ref': comp_ref,
            'name': met.name or met.id,
            'numerical_attributes': {},
            'string_attributes': {},
        })

    # --- modelreactions (skip biomass-like) ----------------------------------
    biomass_rxn_ids = set(BIOMASS_RXNS)
    modelreactions = []
    for r in m.reactions:
        if r.id in biomass_rxn_ids:
            continue
        if r.id.startswith('EX_') or r.id.startswith('DM_') or r.id.startswith('SK_'):
            # Exchange / demand / sink: KBase represents these implicitly via
            # exchanges block (we leave it empty -- the importer typically
            # rebuilds exchanges from media). Include them as model_reactions
            # so they're not lost, but with a transport-like ref.
            pass
        mr = convert_to_kbase_reaction(r, compounds_to_refs)
        if mr is None:
            continue
        # Override modelcompartment_ref to a real compartment
        primary_comp = next(iter({me.compartment for me in r.metabolites}), 'c0')
        mr['modelcompartment_ref'] = compartment_to_ref.get(primary_comp, '~/modelcompartments/id/c0')
        modelreactions.append(mr)

    # --- biomasses ----------------------------------------------------------
    biomasses = []
    for bid, blabel in BIOMASS_RXNS.items():
        if bid not in [r.id for r in m.reactions]:
            continue
        br = m.reactions.get_by_id(bid)
        biomass = {
            'id': bid,
            'name': blabel,
            'other': 0, 'dna': 0, 'rna': 0, 'protein': 0,
            'cellwall': 0, 'lipid': 0, 'cofactor': 0, 'energy': 0,
            'biomasscompounds': [
                {
                    'modelcompound_ref': f'~/modelcompounds/id/{met.id}',
                    'coefficient': float(coef),
                    'edits': {},
                    'gapfill_data': {},
                }
                for met, coef in br.metabolites.items()
            ],
            'edits': {},
            'gapfill_data': {},
            'removedcompounds': [],
            'deleted_compounds': [],
        }
        biomasses.append(biomass)

    # --- top-level kmodel ---------------------------------------------------
    kmodel = {
        'gapfilledcandidates': [],
        'gapgens': [],
        'gapfillings': [],
        'id': object_id,
        'name': object_id,
        # FSP237 genome ref provided by user (2026-06-23): 169876/166/3.
        # Required so the KBase model viewer can resolve gene_NNNN feature
        # references in the GPR display.
        'genome_ref': '169876/166/3',
        'template_ref': '12998/1/2',  # placeholder; safe default for FBAModel template
        'template_refs': ['12998/1/2'],
        'type': 'GenomeScale',
        'source': 'janakakbase: cobra V10 (gapfilled + dedup + dirlock + genes + VLCFA-complete + comp-pruned)',
        'source_id': object_id,
        'biomasses': biomasses,
        'modelcompartments': modelcompartments,
        'modelcompounds': modelcompounds,
        'modelreactions': modelreactions,
    }
    return kmodel


def main():
    print(f'loading model: {MODEL_PATH}')
    m = cobra.io.load_json_model(MODEL_PATH)
    print(f'  {len(m.reactions)} rxns / {len(m.metabolites)} mets / {len(m.genes)} genes')

    print('\nconverting to KBase format (compartments + biomasses fixed) ...')
    kmodel = build_kbase_model(OBJECT_ID, m)
    print(f'  modelreactions   : {len(kmodel["modelreactions"])}')
    print(f'  modelcompounds   : {len(kmodel["modelcompounds"])}')
    print(f'  modelcompartments: {len(kmodel["modelcompartments"])} ({[c["id"] for c in kmodel["modelcompartments"]]})')
    print(f'  biomasses        : {len(kmodel["biomasses"])} ({[b["id"] for b in kmodel["biomasses"]]})')
    print(f'  genome_ref       : {kmodel["genome_ref"]}  (placeholder; re-link in KBase narrative if needed)')
    print(f'  template_ref     : {kmodel["template_ref"]}')

    print('\nconnecting to KBase ...')
    api = cobrakbase.KBaseAPI()
    ws_info = api.get_workspace_info(WS_ID)
    print(f'  workspace        : {WS_ID} = {ws_info[1]} (owned by {ws_info[2]})')

    print(f'\nsaving as {OBJECT_TYPE} -> {WS_ID}/{OBJECT_ID} ...')
    meta = {
        'source':         'COBRA V10 (FSP237 gapfilled + dedup + dirlock + genes + VLCFA-complete + comp-pruned)',
        'n_reactions':    str(len(m.reactions)),
        'n_metabolites':  str(len(m.metabolites)),
        'n_genes':        str(len(m.genes)),
        'github':         'https://github.com/janakagithub/fsp237',
    }
    info = api.save_object(OBJECT_ID, WS_ID, OBJECT_TYPE, kmodel, meta=meta)
    print('\nsaved!')
    # KBaseObjectInfo is an attribute-accessible dataclass-like wrapper
    print(f'  object info: {info}')
    for attr in ('id','name','type','version','ws_id','full_reference'):
        if hasattr(info, attr):
            print(f'    {attr}: {getattr(info, attr)}')


if __name__ == '__main__':
    main()
