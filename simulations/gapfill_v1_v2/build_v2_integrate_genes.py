#!/usr/bin/env python3
"""Phase 4 -- build FSP237 gapfilled Version 2 by integrating gene
assignments onto the Version 1 model.

Inputs:
  - models/fsp237_gapfilled_Version1_noGenes.json     (Phase 1 output)
  - blast/v1_rxn_to_gene_mapping.tsv                  (Phase 2+3 picks)

Logic:
  - For each gapfilled reaction in V1, look up its row in the mapping TSV.
  - If confidence is 'high' or 'medium', assign gene_reaction_rule = fsp237_gene.
  - If confidence is 'low' or 'none', leave GPR blank and record in a flag.
  - Conservative posture: we DO NOT promote a low-confidence match to a GPR
    just to make the model 'complete' -- biology trumps decoration.

For peroxisomal cofactor shuttles (tx_*_xc), the BLAST hits are based on
generic transporter keywords that don't necessarily indicate the *right*
membrane transporter. We treat them as 'unassigned' in V2 by default;
flip the FLAG_SHUTTLES_UNASSIGNED switch if you want to use the BLAST picks
verbatim instead.

Saves:
  - models/fsp237_gapfilled_Version2_gapfill_genes_integrated.json
  - reports/v2_gpr_assignments.tsv  (audit trail: rxn_id, gpr, confidence, note)
"""
import csv
import os

import cobra

BASE     = '/home/janakae/fungalTemplate/imm904CobraModel'
HERE     = f'{BASE}/simulations/gapfill_v1_v2'
V1_MODEL = f'{HERE}/models/fsp237_gapfilled_Version1_noGenes.json'
V2_MODEL = f'{HERE}/models/fsp237_gapfilled_Version2_gapfill_genes_integrated.json'
MAP_TSV  = f'{HERE}/blast/v1_rxn_to_gene_mapping.tsv'
AUDIT    = f'{HERE}/reports/v2_gpr_assignments.tsv'

# Conservative: do NOT integrate shuttle gene picks (low biological specificity)
FLAG_SHUTTLES_UNASSIGNED = True


def main():
    print(f'loading: {V1_MODEL}')
    m = cobra.io.load_json_model(V1_MODEL)
    print(f'  {len(m.reactions)} rxns / {len(m.metabolites)} mets / {len(m.genes)} genes')

    # Load mapping
    picks = {}
    with open(MAP_TSV) as fh:
        rdr = csv.DictReader(fh, delimiter='\t')
        for row in rdr:
            picks[row['rxn_id']] = row

    print(f'\nloaded {len(picks)} pick rows')

    # Apply GPRs
    audit = []
    n_assigned = 0
    n_skipped  = 0
    n_low      = 0
    for rxn_id, p in picks.items():
        is_shuttle = rxn_id.startswith('tx_')
        gene = p['fsp237_gene']
        conf = p['confidence']
        decision = 'assign'
        reason = ''
        if conf in ('none',) or not gene:
            decision = 'skip-no-match'
            reason = 'No confident BLAST hit; leave GPR empty for manual review'
        elif conf == 'low':
            decision = 'skip-low-conf'
            reason = f'BLAST pident={p["pident"]} qcov={p["qcovs"]} below medium threshold; flag for review'
        elif is_shuttle and FLAG_SHUTTLES_UNASSIGNED:
            decision = 'skip-shuttle'
            reason = ('Generic transporter keywords -- BLAST hit may not be the '
                       'correct peroxisomal carrier. Left unassigned for manual review.')
        if decision == 'assign':
            try:
                r = m.reactions.get_by_id(rxn_id)
                r.gene_reaction_rule = gene
                n_assigned += 1
            except KeyError:
                decision = 'skip-rxn-not-in-model'
                reason = f'Reaction {rxn_id} not present in V1 model'
                n_skipped += 1
        else:
            if conf == 'low': n_low += 1
            else:             n_skipped += 1

        audit.append({
            'rxn_id': rxn_id,
            'decision': decision,
            'gpr_assigned': gene if decision == 'assign' else '',
            'confidence': conf,
            'source_locus_tag': p['source_locus_tag'],
            'source_product': p['source_product'],
            'source_strain': p['source_strain'],
            'pident': p['pident'], 'qcovs': p['qcovs'],
            'evalue': p['evalue'], 'bitscore': p['bitscore'],
            'reason': reason or 'high-confidence BLAST hit; assigned as GPR',
        })

    print(f'\nassigned : {n_assigned}')
    print(f'skipped  : {n_skipped}  (no match / shuttle / not in model)')
    print(f'low-conf : {n_low}')

    cobra.io.save_json_model(m, V2_MODEL)
    print(f'\nsaved: {V2_MODEL}')

    with open(AUDIT, 'w', newline='') as fh:
        cols = ['rxn_id','decision','gpr_assigned','confidence',
                'source_locus_tag','source_product','source_strain',
                'pident','qcovs','evalue','bitscore','reason']
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        w.writeheader()
        for row in audit: w.writerow(row)
    print(f'saved: {AUDIT}')


if __name__ == '__main__':
    main()
