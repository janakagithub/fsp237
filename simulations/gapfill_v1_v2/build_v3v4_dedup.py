#!/usr/bin/env python3
"""V3 / V4 -- model-wide exact-duplicate cleanup.

Inputs:
  V1: simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version1_noGenes.json
  V2: simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version2_gapfill_genes_integrated.json

For each input we:
  1. Find every group of reactions that share the same compartment-aware
     stoichiometry (treating reverse-written variants as equivalent if both
     are reversible or if either spans both directions).
  2. Pick one canonical representative per group using a precedence rule:
       a. Prefer a 'rxn'-prefixed ModelSEED id over BiGG-style aliases
          (ALCD2ir_c0, HDCAt_c0, etc.).
       b. Among ModelSEED ids, prefer the lower numeric suffix (rxn00102
          beats rxn30392) -- matches the dedup policy from
          BuildMinimalFSP237Model.ipynb.
       c. Prefer 'rxn' over 'frxn' (canonical over fix).
       d. Tiebreak alphabetically.
  3. Union-merge the GPRs of the duplicates onto the kept rxn (collected
     as ' or '-joined unique gene tokens). Bounds are widened to the
     superset of the duplicates' bounds.
  4. Remove the duplicates.

Outputs:
  models/fsp237_gapfilled_Version3_dedup_noGenes.json          (from V1)
  models/fsp237_gapfilled_Version4_dedup_genes_integrated.json (from V2)
  reports/v3_dedup_log.tsv  -- one row per dedup decision, with kept/dropped
                                ids, transferred GPR tokens, and rationale.
"""
import csv
import os
import re
from collections import defaultdict

import cobra

BASE = '/home/janakae/fungalTemplate/imm904CobraModel'
HERE = f'{BASE}/simulations/gapfill_v1_v2'
V1   = f'{HERE}/models/fsp237_gapfilled_Version1_noGenes.json'
V2   = f'{HERE}/models/fsp237_gapfilled_Version2_gapfill_genes_integrated.json'
V3   = f'{HERE}/models/fsp237_gapfilled_Version3_dedup_noGenes.json'
V4   = f'{HERE}/models/fsp237_gapfilled_Version4_dedup_genes_integrated.json'
LOG_V3 = f'{HERE}/reports/v3_dedup_log.tsv'
LOG_V4 = f'{HERE}/reports/v4_dedup_log.tsv'


# Reactions to PROTECT from dedup -- functional duplicates we deliberately
# keep as separate objectives / role markers. Each PROTECTED rxn is removed
# from its dedup group before precedence is applied, so both the protected
# rxn and the regular canonical rxn survive.
#
# `bio1` is the ATP-yield objective used by extend_biomass.py,
# build_atp_safe_site.py, build_escher_maps.py, etc. -- stoichiometrically
# identical to rxn00062_c0 (NGAM) but kept as a separate reaction so the
# 30/2 ATP yield test is still pointable. Per user instruction (2026-06-23).
PROTECTED_IDS = {
    'bio1',
}

MSEED_RE = re.compile(r'^rxn(\d+)(.*)$')   # rxn00102_c0 etc.

def _is_modelseed(rid):
    return rid.startswith('rxn') and not rid.startswith('rxnL_')

def _modelseed_num(rid):
    m = MSEED_RE.match(rid)
    return int(m.group(1)) if m else None

def pick_canonical(ids):
    """Apply precedence rules; return the chosen rxn_id."""
    # 1. Prefer ModelSEED-prefixed
    mseed = [r for r in ids if _is_modelseed(r)]
    if mseed:
        # 2. Lowest numeric suffix
        mseed.sort(key=lambda r: (_modelseed_num(r) is None,
                                    _modelseed_num(r) or 1_000_000_000, r))
        return mseed[0]
    # 3. Prefer 'rxn' over 'frxn'
    frxn = [r for r in ids if r.startswith('rxn')]
    if frxn:
        return sorted(frxn)[0]
    # 4. Alphabetical fallback
    return sorted(ids)[0]


def _stoich_sig(rxn):
    """Compartment-aware stoichiometry signature (forward + reverse forms)."""
    items = tuple(sorted((m.id, c) for m, c in rxn.metabolites.items()))
    rev   = tuple(sorted((mid, -c) for mid, c in items))
    return frozenset([items, rev])


def merge_gprs(kept_gpr, dropped_gpr):
    """Union the 'or'-style GPRs, deduping gene tokens, preserving order."""
    def tokens(s):
        if not s: return []
        return [t.strip() for t in re.split(r'\s+or\s+', s) if t.strip()]
    seen = set()
    out  = []
    for t in tokens(kept_gpr) + tokens(dropped_gpr):
        if t not in seen:
            seen.add(t); out.append(t)
    return ' or '.join(out)


def widen_bounds(kept, dropped):
    lb = min(kept.lower_bound, dropped.lower_bound)
    ub = max(kept.upper_bound, dropped.upper_bound)
    return lb, ub


def dedup_model(in_path, out_path, log_path):
    print(f'\n========== DEDUP {os.path.basename(in_path)} ==========')
    m = cobra.io.load_json_model(in_path)
    print(f'  start: {len(m.reactions)} rxns, {len(m.metabolites)} mets, {len(m.genes)} genes')

    # Group reactions by signature (compartment-aware)
    sig_to_ids = defaultdict(list)
    for r in m.reactions:
        for s in _stoich_sig(r):
            sig_to_ids[s].append(r.id)

    # Build groups (union-find via shared sigs)
    visited = set()
    groups = []
    for r in m.reactions:
        if r.id in visited: continue
        members = set()
        for s in _stoich_sig(r):
            for other_id in sig_to_ids[s]:
                # require same compartment AS_A_SET so we don't merge e.g. c0+m0
                other = m.reactions.get_by_id(other_id)
                if {x.compartment for x in other.metabolites} == {x.compartment for x in r.metabolites}:
                    members.add(other_id)
        if len(members) > 1:
            # Drop any PROTECTED_IDS from the dedup group so they survive.
            members -= PROTECTED_IDS
        if len(members) > 1:
            groups.append(sorted(members))
            visited.update(members)

    print(f'  found {len(groups)} duplicate groups; {sum(len(g) for g in groups)} reactions involved')

    log = []
    n_dropped = 0
    for g in groups:
        keeper = pick_canonical(g)
        keeper_rxn = m.reactions.get_by_id(keeper)
        dropped = [rid for rid in g if rid != keeper]
        kept_gpr_old = keeper_rxn.gene_reaction_rule or ''
        new_gpr = kept_gpr_old
        new_bounds = (keeper_rxn.lower_bound, keeper_rxn.upper_bound)
        for drop_id in dropped:
            drop_rxn = m.reactions.get_by_id(drop_id)
            new_gpr = merge_gprs(new_gpr, drop_rxn.gene_reaction_rule or '')
            new_bounds = widen_bounds(
                type('B', (), {'lower_bound': new_bounds[0], 'upper_bound': new_bounds[1]}),
                drop_rxn,
            )
        # Apply merged GPR + widened bounds, then remove duplicates
        gpr_changed = new_gpr != kept_gpr_old
        bounds_changed = new_bounds != (keeper_rxn.lower_bound, keeper_rxn.upper_bound)
        if gpr_changed:
            keeper_rxn.gene_reaction_rule = new_gpr
        if bounds_changed:
            keeper_rxn.lower_bound, keeper_rxn.upper_bound = new_bounds
        log.append({
            'group_size': len(g),
            'kept': keeper,
            'dropped': ';'.join(dropped),
            'kept_name': keeper_rxn.name or '',
            'kept_gpr_before': kept_gpr_old,
            'kept_gpr_after': new_gpr,
            'gpr_changed': str(gpr_changed),
            'bounds_changed': str(bounds_changed),
            'kept_bounds': f'{keeper_rxn.lower_bound},{keeper_rxn.upper_bound}',
            'reason': 'ModelSEED-precedence + lowest-numeric-suffix dedup',
        })
        m.remove_reactions(dropped, remove_orphans=False)
        n_dropped += len(dropped)

    print(f'  removed {n_dropped} duplicates')
    print(f'  end  : {len(m.reactions)} rxns, {len(m.metabolites)} mets, {len(m.genes)} genes')

    cobra.io.save_json_model(m, out_path)
    print(f'  saved: {out_path}')

    with open(log_path, 'w', newline='') as fh:
        cols = ['group_size','kept','dropped','kept_name','kept_gpr_before',
                'kept_gpr_after','gpr_changed','bounds_changed','kept_bounds','reason']
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        w.writeheader()
        for row in log: w.writerow(row)
    print(f'  log  : {log_path}')

    return n_dropped


def main():
    n1 = dedup_model(V1, V3, LOG_V3)
    n2 = dedup_model(V2, V4, LOG_V4)
    print(f'\nSUMMARY: V1 -> V3 dropped {n1}; V2 -> V4 dropped {n2}')


if __name__ == '__main__':
    main()
