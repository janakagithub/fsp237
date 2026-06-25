#!/usr/bin/env python3
"""Reaction-deduplication pass extracted from BuildMinimalFSP237Model.ipynb
(cell 15).

This is the FIRST dedup pass in the FSP237 build chain. It is exact-
stoichiometry (signs matter, so reverse-written duplicates are NOT
collapsed) and uses a simple "lowest-numeric-suffix ModelSEED ID wins"
precedence rule.

When to use this vs the V3/V4 dedup:
  - This script ── runs during/after the initial KBase + Excel + iMM904
                    merge, when you have many same-direction duplicates
                    from different upstream sources. Conservative: only
                    collapses reactions whose stoichiometry hash matches
                    EXACTLY (sign-sensitive).
  - simulations/gapfill_v1_v2/build_v3v4_dedup.py
                 ── runs on the gap-fill output. Handles BOTH the forward
                    and reverse stoichiometry hashes (so a reaction written
                    "A->B" collapses with one written "B->A" if they are
                    biologically the same). Also has a PROTECTED_IDS set
                    so `bio1` is preserved as a separate ATP-yield objective.
                    Use this for the final cleanup.

Both scripts:
  - Skip H+ when building the stoichiometry hash (cofactor balancing
    artifacts shouldn't separate otherwise-identical reactions).
  - Prefer ModelSEED `rxn*` IDs over BiGG/yeast aliases.
  - Among ModelSEED IDs, the lower numeric suffix wins (rxn00102 beats
    rxn30392).
  - Union-merge GPRs from collapsed duplicates onto the kept reaction.
  - Widen bounds to the superset of the duplicates' bounds.
  - Remove orphan metabolites left behind.

Usage:
  python dedup_initial_build.py <input_model.json> [<output_model.json>]
  # default output: <input>_dedup.json next to the input

  # also accepts --log <path> to write the per-decision audit TSV.
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

import cobra


# Cofactors to skip when hashing stoichiometry (so a charge-balanced
# H+ difference doesn't keep otherwise-identical reactions apart).
COFACTORS_BASE = {'cpd00067'}   # H+

MS_RE = re.compile(r'^rxn(\d+)')


def canon_stoich(rxn):
    """Canonical (sign-sensitive) stoichiometry hash, excluding H+.

    Returns a tuple of (metabolite_id, signed_coefficient) sorted by id.
    Reverse-written variants (negated signs) produce DIFFERENT hashes by
    design -- this script is conservative.
    """
    items = []
    for met, coef in rxn.metabolites.items():
        if met.id.split('_')[0] in COFACTORS_BASE:
            continue
        items.append((met.id, coef))
    return tuple(sorted(items))


def ms_num(rid):
    """Return the integer suffix of a ModelSEED-style id, or None."""
    m = MS_RE.match(rid)
    return int(m.group(1)) if m else None


def pick_keeper(ids):
    """Lower-numeric-suffix ModelSEED ID wins; non-MS IDs are last."""
    ms_ids = sorted([r for r in ids if ms_num(r) is not None], key=ms_num)
    non_ms = [r for r in ids if ms_num(r) is None]
    if ms_ids:
        return ms_ids[0], ms_ids[1:] + non_ms
    return ids[0], ids[1:]


def dedup(model):
    """Run the dedup in-place on `model`. Returns (log_rows, n_orphans_removed)."""
    dup_groups = defaultdict(list)
    for r in model.reactions:
        k = canon_stoich(r)
        if not k:
            continue
        dup_groups[k].append(r.id)
    dup_groups = {k: v for k, v in dup_groups.items() if len(v) > 1}
    print(f'  duplicate groups found: {len(dup_groups)}')

    log = []
    for ids in dup_groups.values():
        keep, remove = pick_keeper(ids)
        kp = model.reactions.get_by_id(keep)
        for rid in remove:
            rm = model.reactions.get_by_id(rid)

            # Conservative GPR merge: only adopt the removed reaction's GPR
            # if the kept reaction has none. (build_v3v4_dedup.py does a
            # true union-merge; here we mirror the original notebook logic.)
            merged_gpr = False
            if not kp.gene_reaction_rule and rm.gene_reaction_rule:
                kp.gene_reaction_rule = rm.gene_reaction_rule
                merged_gpr = True

            new_lb = min(kp.lower_bound, rm.lower_bound)
            new_ub = max(kp.upper_bound, rm.upper_bound)
            widened = (new_lb, new_ub) != (kp.lower_bound, kp.upper_bound)
            if widened:
                kp.lower_bound, kp.upper_bound = new_lb, new_ub

            log.append({
                'removed_rxn': rid,
                'kept_rxn': keep,
                'removed_name': rm.name,
                'kept_name': kp.name,
                'removed_eq': rm.build_reaction_string(),
                'kept_eq': kp.build_reaction_string(),
                'removed_lb': rm.lower_bound,
                'removed_ub': rm.upper_bound,
                'kept_lb': kp.lower_bound,
                'kept_ub': kp.upper_bound,
                'removed_genes': ';'.join(g.id for g in rm.genes),
                'kept_genes': ';'.join(g.id for g in kp.genes),
                'merged_gpr_from_removed': merged_gpr,
                'widened_bounds_to': f'[{new_lb},{new_ub}]' if widened else '',
            })

    model.remove_reactions([e['removed_rxn'] for e in log])
    orphans = [m for m in model.metabolites if len(m.reactions) == 0]
    n_orphans = len(orphans)
    if orphans:
        model.remove_metabolites(orphans)
    return log, n_orphans


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input_model', help='Path to COBRA JSON model to dedup')
    ap.add_argument('output_model', nargs='?', default=None,
                    help='Output path (default: <input>_dedup.json)')
    ap.add_argument('--log', default=None,
                    help='Per-decision audit TSV (default: <input>_dedup_log.tsv)')
    args = ap.parse_args()

    in_path = args.input_model
    out_path = args.output_model or in_path.replace('.json', '_dedup.json')
    log_path = args.log or in_path.replace('.json', '_dedup_log.tsv')

    print(f'loading: {in_path}')
    m = cobra.io.load_json_model(in_path)
    print(f'  start: {len(m.reactions)} rxns, {len(m.metabolites)} mets, {len(m.genes)} genes')

    log, n_orph = dedup(m)
    print(f'  removed {len(log)} duplicate reactions, {n_orph} orphan metabolites')
    print(f'  end  : {len(m.reactions)} rxns, {len(m.metabolites)} mets, {len(m.genes)} genes')

    # Smoke test
    sol = m.optimize()
    print(f'  biomass after dedup: {sol.objective_value or 0:.6g}  (status={sol.status})')

    cobra.io.save_json_model(m, out_path)
    print(f'  saved: {out_path}')

    if log:
        with open(log_path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(log[0].keys()), delimiter='\t')
            w.writeheader()
            for row in log:
                w.writerow(row)
        print(f'  log  : {log_path} ({len(log)} entries)')


if __name__ == '__main__':
    main()
