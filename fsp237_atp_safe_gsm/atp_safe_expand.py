"""Helpers for ATP-safe expansion of the FSP237 CMM into a GSM.

The CMM (KBase 169876/170/5) is curated to give exactly 30 ATP/glucose aerobic
and 2 ATP/glucose anaerobic. We expand it by adding reactions from the larger
fsp237 GSM, rejecting (blacklisting) any reaction that pushes ATP yield above
the curated value -- that's the signature of an introduced futile cycle.

Workflow:
    1. minimal_media(model, glc=-1) sets the standard ATP-yield test bounds.
    2. atp_yield(model) returns (aerobic_yield, anaerobic_yield).
    3. expand_atp_safe(cmm, gsm, ...) iteratively adds the GSM-unique reactions
       in batches, bisecting on failure, and returns (safe_added, blacklisted).
    4. test_biomass(model, biomass_rxn) returns the biomass FBA value.
    5. minimum_blacklist_for_biomass(...) re-adds the smallest blacklist subset
       that lets biomass grow above a threshold.
"""
from __future__ import annotations

import copy
import time
from dataclasses import dataclass, field
from typing import Iterable

import cobra


# ---- media presets -----------------------------------------------------------

# These exchanges stay open in BOTH directions in every test (mass-balance
# necessities + standard fungal minimal-media nutrients). Verified that adding
# the trace ions / sulfate / NH3 does NOT change the CMM's ATP yield (still
# 30 aerobic / 2 anaerobic), so the wider set is safe to use throughout.
# Everything else gets closed to uptake; only what we open below can enter.
ALWAYS_OPEN = {
    'cpd00001_e0',  # H2O
    'cpd00009_e0',  # Phosphate (Pi)
    'cpd00011_e0',  # CO2
    'cpd00067_e0',  # H+
    'cpd00013_e0',  # NH3
    'cpd00048_e0',  # Sulfate (needed for biomass: Cys, Met)
    'cpd10515_e0',  # Fe2+
    'cpd00971_e0',  # Na+
    'cpd00205_e0',  # K+
}

GLC_EX = 'EX_cpd00027_e0'
O2_EX  = 'EX_cpd00007_e0'


def minimal_media(model: cobra.Model, glc_uptake: float = 1.0,
                  aerobic: bool = True, o2_uptake: float = 1000.0) -> None:
    """In-place: close every exchange, then open the standard test set.

    glc_uptake: positive number; will be applied as lower_bound = -glc_uptake.
    aerobic: if True, O2 is opened to -o2_uptake; else O2 stays closed.
    """
    for ex in model.exchanges:
        met_id = next(iter(ex.metabolites.keys())).id
        if met_id in ALWAYS_OPEN:
            ex.lower_bound = -1000; ex.upper_bound = 1000
        else:
            ex.lower_bound = 0; ex.upper_bound = 1000

    glc = model.reactions.get_by_id(GLC_EX)
    glc.lower_bound = -glc_uptake; glc.upper_bound = 0

    o2 = model.reactions.get_by_id(O2_EX)
    if aerobic:
        o2.lower_bound = -o2_uptake; o2.upper_bound = 0
    else:
        o2.lower_bound = 0; o2.upper_bound = 0


def atp_yield(model: cobra.Model, atp_rxn: str, glc_uptake: float = 1.0,
              tol: float = 1e-6) -> tuple[float, float]:
    """Return (aerobic, anaerobic) ATP yield per glucose, given the model already
    has the ATP-demand reaction set as objective via `atp_rxn`.

    Uses a context so the media reset doesn't leak.
    """
    with model:
        # Force objective onto the ATP-demand reaction
        for r in model.reactions:
            r.objective_coefficient = 1.0 if r.id == atp_rxn else 0.0
        minimal_media(model, glc_uptake=glc_uptake, aerobic=True)
        sol_a = model.optimize()
        a = sol_a.objective_value if sol_a.objective_value is not None else float('nan')

        minimal_media(model, glc_uptake=glc_uptake, aerobic=False)
        sol_n = model.optimize()
        n = sol_n.objective_value if sol_n.objective_value is not None else float('nan')

    return a / glc_uptake, n / glc_uptake


# ---- candidate selection -----------------------------------------------------

def gsm_unique_reactions(cmm: cobra.Model, gsm: cobra.Model) -> list[str]:
    """Reactions in GSM not present in CMM (by reaction id)."""
    cmm_ids = {r.id for r in cmm.reactions}
    return [r.id for r in gsm.reactions if r.id not in cmm_ids]


# ---- shared-reaction handling ------------------------------------------------

def _swap_to_gsm(model: cobra.Model, gsm: cobra.Model, rid: str) -> tuple:
    """Replace model's reaction `rid` with the GSM's version (stoich + bounds).
    Returns (old_bounds, old_stoich) so the caller can revert."""
    cr = model.reactions.get_by_id(rid)
    gr = gsm.reactions.get_by_id(rid)
    old_bounds = cr.bounds
    old_stoich = dict(cr.metabolites)
    cr.subtract_metabolites(cr.metabolites)
    missing = [m for m in gr.metabolites if m.id not in [x.id for x in model.metabolites]]
    if missing:
        model.add_metabolites([cobra.Metabolite(id=m.id, name=m.name, formula=m.formula,
                                                  charge=m.charge, compartment=m.compartment)
                                for m in missing])
    cr.add_metabolites({model.metabolites.get_by_id(m.id): c for m, c in gr.metabolites.items()})
    cr.bounds = gr.bounds
    return old_bounds, old_stoich


def _revert_swap(model: cobra.Model, rid: str, old_bounds, old_stoich):
    cr = model.reactions.get_by_id(rid)
    cr.subtract_metabolites(cr.metabolites)
    cr.add_metabolites(old_stoich)
    cr.bounds = old_bounds


def find_mismatched_shared(cmm: cobra.Model, gsm: cobra.Model) -> list[str]:
    """Return ids of reactions present in both models but with different
    bounds and/or stoichiometry."""
    out = []
    cmm_rxns = {r.id: r for r in cmm.reactions}
    for rid, gr in {r.id: r for r in gsm.reactions}.items():
        if rid not in cmm_rxns: continue
        cr = cmm_rxns[rid]
        cs = {m.id: c for m, c in cr.metabolites.items()}
        gs = {m.id: c for m, c in gr.metabolites.items()}
        if cs != gs or cr.bounds != gr.bounds:
            out.append(rid)
    return out


def minimum_shared_swaps_for_biomass(cmm: cobra.Model, gsm: cobra.Model,
                                      biomass_rxn: str, atp_rxn: str,
                                      baseline_yields: tuple[float, float],
                                      glc_uptake_atp: float = 1.0,
                                      glc_uptake_bio: float = 5.0,
                                      biomass_threshold: float = 1e-6,
                                      atp_tol: float = 1e-3,
                                      use_gsm_media: bool = True,
                                      verbose: bool = True) -> tuple[list[str], list[str]]:
    """CMM has precedence. Only swap shared reactions to GSM version when
    biomass cannot grow otherwise. Returns (swapped, kept_cmm).

    Algorithm:
        1. Test biomass with NO swaps. If > threshold, return ([], [all_mismatched]).
        2. Try swapping ALL mismatched shared rxns. If biomass still 0, no swap
           combination helps; return ([], [all_mismatched]).
        3. Otherwise bisect to find the smallest subset that enables biomass.
        4. Among the candidate swap set, also verify ATP yield is preserved;
           if a swap breaks ATP yield, keep it as CMM.

    Uses GSM-style media for biomass test if use_gsm_media=True.
    """
    mismatched = find_mismatched_shared(cmm, gsm)
    if verbose:
        print(f'shared rxns with bounds/stoich mismatch: {len(mismatched)}')

    def _set_media_for_biomass(model):
        if use_gsm_media:
            # copy GSM exchange bounds onto target
            gsm_bounds = {r.id: (r.lower_bound, r.upper_bound) for r in gsm.exchanges}
            for ex in model.exchanges:
                if ex.id in gsm_bounds:
                    ex.lower_bound, ex.upper_bound = gsm_bounds[ex.id]
                else:
                    ex.lower_bound = 0
        else:
            minimal_media(model, glc_uptake=glc_uptake_bio, aerobic=True)

    def _test_bm(model):
        with model:
            _set_media_for_biomass(model)
            for r in model.reactions:
                r.objective_coefficient = 1.0 if r.id == biomass_rxn else 0.0
            sol = model.optimize()
            return sol.objective_value if sol.objective_value is not None else 0.0

    # Step 1: try without any swap
    bm0 = _test_bm(cmm)
    if bm0 >= biomass_threshold:
        if verbose:
            print(f'biomass grows without any shared-rxn swap: {bm0:.6f}')
        return [], mismatched

    # Step 2: swap all, check if it grows
    saves = []  # (rid, old_bounds, old_stoich) tuples for revert
    for rid in mismatched:
        ob, os_ = _swap_to_gsm(cmm, gsm, rid)
        saves.append((rid, ob, os_))
    bm_all = _test_bm(cmm)
    if verbose:
        print(f'biomass with ALL shared rxns swapped: {bm_all:.6f}')
    if bm_all < biomass_threshold:
        # No swap combination helps; revert all
        for rid, ob, os_ in reversed(saves):
            _revert_swap(cmm, rid, ob, os_)
        return [], mismatched

    # Step 3: greedy reduction -- revert one swap at a time; keep if biomass holds
    if verbose:
        print(f'greedily reducing swaps to minimum subset...')
    needed = list(mismatched)
    for rid in list(mismatched):
        # try reverting this one
        cr = cmm.reactions.get_by_id(rid)
        cur_bounds = cr.bounds
        cur_stoich = dict(cr.metabolites)
        # find the original
        orig = next((ob, os_) for r, ob, os_ in saves if r == rid)
        _revert_swap(cmm, rid, orig[0], orig[1])
        bm = _test_bm(cmm)
        if bm >= biomass_threshold:
            needed.remove(rid)  # don't need this swap
        else:
            # restore swap
            _swap_to_gsm(cmm, gsm, rid)
    if verbose:
        print(f'minimum swap set for biomass: {len(needed)} -> {needed}')

    # Step 4: verify ATP yield is preserved with the minimum-needed set
    yields = atp_yield(cmm, atp_rxn, glc_uptake=glc_uptake_atp)
    if verbose:
        print(f'ATP yield with minimum swap set: {yields[0]:.3f}/{yields[1]:.3f}')
    if (abs(yields[0] - baseline_yields[0]) > atp_tol or
        abs(yields[1] - baseline_yields[1]) > atp_tol):
        if verbose:
            print(f'WARNING: swap set breaks ATP yield; this is a tradeoff -- '
                  f'biomass requires reactions that perturb ATP yield')
    kept_cmm = [rid for rid in mismatched if rid not in needed]
    return needed, kept_cmm


# ---- core: add a reaction (and its missing metabolites) into a target model --

def add_reactions_from(target: cobra.Model, source: cobra.Model,
                       rxn_ids: Iterable[str]) -> list[str]:
    """Copy reactions from source into target, adding any missing metabolites
    along the way. Returns the list of ids actually added (skips already-present)."""
    added = []
    tgt_rxn_ids = {r.id for r in target.reactions}
    tgt_met_ids = {m.id for m in target.metabolites}

    # First pass: collect & add missing metabolites
    new_mets = {}
    for rid in rxn_ids:
        if rid in tgt_rxn_ids:
            continue
        src_r = source.reactions.get_by_id(rid)
        for m in src_r.metabolites:
            if m.id not in tgt_met_ids and m.id not in new_mets:
                new_m = cobra.Metabolite(
                    id=m.id, name=m.name, formula=m.formula,
                    charge=m.charge, compartment=m.compartment,
                )
                new_m.annotation = dict(m.annotation)
                new_mets[m.id] = new_m
    if new_mets:
        target.add_metabolites(list(new_mets.values()))
        tgt_met_ids.update(new_mets)

    # Second pass: add reactions (mapping their metabolites onto target's)
    new_rxns = []
    for rid in rxn_ids:
        if rid in tgt_rxn_ids:
            continue
        src_r = source.reactions.get_by_id(rid)
        new_r = cobra.Reaction(id=rid, name=src_r.name,
                                lower_bound=src_r.lower_bound,
                                upper_bound=src_r.upper_bound)
        new_r.gene_reaction_rule = src_r.gene_reaction_rule
        new_r.subsystem = src_r.subsystem
        new_r.annotation = dict(src_r.annotation)
        # Map metabolites to the target's instances
        stoich = {target.metabolites.get_by_id(m.id): coef
                  for m, coef in src_r.metabolites.items()}
        new_r.add_metabolites(stoich)
        new_rxns.append(new_r)
        added.append(rid)
    if new_rxns:
        target.add_reactions(new_rxns)
    return added


def remove_reactions(model: cobra.Model, rxn_ids: Iterable[str]) -> None:
    rxns = [model.reactions.get_by_id(rid) for rid in rxn_ids
            if rid in [r.id for r in model.reactions]]
    if rxns:
        model.remove_reactions(rxns, remove_orphans=False)


# ---- iterative ATP-safe expansion -------------------------------------------

@dataclass
class ExpansionResult:
    safe_added: list[str] = field(default_factory=list)
    blacklisted: list[str] = field(default_factory=list)
    baseline_yields: tuple[float, float] = (0.0, 0.0)
    final_yields:    tuple[float, float] = (0.0, 0.0)
    batch_log: list[dict] = field(default_factory=list)


def _yields_ok(actual: tuple[float, float], baseline: tuple[float, float],
               tol: float) -> bool:
    """ATP yield must not increase above baseline (only adding rxns; can't decrease).
    Allow a small numerical slack."""
    a_ok = (actual[0] is not None and not _isnan(actual[0])
            and actual[0] <= baseline[0] + tol)
    n_ok = (actual[1] is not None and not _isnan(actual[1])
            and actual[1] <= baseline[1] + tol)
    # also reject NaNs / infeasibles
    if _isnan(actual[0]) or _isnan(actual[1]):
        return False
    return a_ok and n_ok


def _isnan(x):
    try:
        return x != x
    except Exception:
        return True


def expand_atp_safe(cmm: cobra.Model, gsm: cobra.Model, atp_rxn: str,
                    candidates: list[str] | None = None,
                    batch_size: int = 50, tol: float = 1e-3,
                    glc_uptake: float = 1.0,
                    verbose: bool = True) -> ExpansionResult:
    """Iteratively add GSM-unique reactions to CMM. After each batch, recheck
    the (aerobic, anaerobic) ATP yields. If they exceed the baseline + tol,
    bisect the batch to find the offenders, blacklist them, and continue.

    Modifies cmm in place.
    """
    result = ExpansionResult()
    result.baseline_yields = atp_yield(cmm, atp_rxn, glc_uptake=glc_uptake)
    if verbose:
        print(f'baseline ATP yields (aerobic, anaerobic): '
              f'{result.baseline_yields[0]:.3f}, {result.baseline_yields[1]:.3f}')

    if candidates is None:
        candidates = gsm_unique_reactions(cmm, gsm)
    if verbose:
        print(f'candidate reactions to consider: {len(candidates)}')

    # Slice into batches
    remaining = list(candidates)
    batch_num = 0
    t0 = time.time()
    while remaining:
        batch_num += 1
        batch = remaining[:batch_size]
        remaining = remaining[batch_size:]
        safe_in_batch, bad_in_batch = _try_batch(
            cmm, gsm, batch, atp_rxn, glc_uptake, result.baseline_yields, tol)
        result.safe_added.extend(safe_in_batch)
        result.blacklisted.extend(bad_in_batch)
        elapsed = time.time() - t0
        result.batch_log.append({
            'batch': batch_num, 'size': len(batch),
            'safe': len(safe_in_batch), 'bad': len(bad_in_batch),
            'cumulative_safe': len(result.safe_added),
            'cumulative_bad':  len(result.blacklisted),
            'elapsed_s': round(elapsed, 1),
        })
        if verbose:
            print(f'  batch {batch_num:3d}: tried {len(batch):3d}, '
                  f'safe={len(safe_in_batch):3d}, '
                  f'bad={len(bad_in_batch):3d} | '
                  f'cum safe={len(result.safe_added)}, '
                  f'cum bad={len(result.blacklisted)} | '
                  f'{elapsed:.1f}s elapsed')

    result.final_yields = atp_yield(cmm, atp_rxn, glc_uptake=glc_uptake)
    if verbose:
        print(f'final ATP yields (aerobic, anaerobic): '
              f'{result.final_yields[0]:.3f}, {result.final_yields[1]:.3f}')
    return result


def _try_batch(cmm, gsm, batch, atp_rxn, glc_uptake, baseline, tol):
    """Try adding a whole batch. If OK, great. If not, bisect to find offenders."""
    add_reactions_from(cmm, gsm, batch)
    yields = atp_yield(cmm, atp_rxn, glc_uptake=glc_uptake)
    if _yields_ok(yields, baseline, tol):
        return list(batch), []

    # Roll back the entire batch, then bisect
    remove_reactions(cmm, batch)
    safe, bad = [], []
    _bisect(cmm, gsm, batch, atp_rxn, glc_uptake, baseline, tol, safe, bad)
    return safe, bad


def _bisect(cmm, gsm, batch, atp_rxn, glc_uptake, baseline, tol, safe, bad):
    """Recursive bisection. Adds whatever can be safely added; blacklists what can't."""
    if not batch:
        return
    if len(batch) == 1:
        rid = batch[0]
        add_reactions_from(cmm, gsm, [rid])
        yields = atp_yield(cmm, atp_rxn, glc_uptake=glc_uptake)
        if _yields_ok(yields, baseline, tol):
            safe.append(rid)
        else:
            remove_reactions(cmm, [rid])
            bad.append(rid)
        return

    # Try first half
    half = len(batch) // 2
    left, right = batch[:half], batch[half:]
    add_reactions_from(cmm, gsm, left)
    yields = atp_yield(cmm, atp_rxn, glc_uptake=glc_uptake)
    if _yields_ok(yields, baseline, tol):
        safe.extend(left)
    else:
        remove_reactions(cmm, left)
        _bisect(cmm, gsm, left, atp_rxn, glc_uptake, baseline, tol, safe, bad)
    add_reactions_from(cmm, gsm, right)
    yields = atp_yield(cmm, atp_rxn, glc_uptake=glc_uptake)
    if _yields_ok(yields, baseline, tol):
        safe.extend(right)
    else:
        remove_reactions(cmm, right)
        _bisect(cmm, gsm, right, atp_rxn, glc_uptake, baseline, tol, safe, bad)


# ---- biomass test + minimum blacklist re-introduction ------------------------

def test_biomass(model: cobra.Model, biomass_rxn: str,
                 glc_uptake: float = 1.0, aerobic: bool = True) -> float:
    with model:
        for r in model.reactions:
            r.objective_coefficient = 1.0 if r.id == biomass_rxn else 0.0
        minimal_media(model, glc_uptake=glc_uptake, aerobic=aerobic)
        sol = model.optimize()
        return sol.objective_value if sol.objective_value is not None else 0.0


def minimum_blacklist_for_biomass(model: cobra.Model, gsm: cobra.Model,
                                   blacklist: list[str], biomass_rxn: str,
                                   atp_rxn: str, baseline_yields: tuple[float, float],
                                   glc_uptake: float = 1.0,
                                   biomass_threshold: float = 1e-6,
                                   verbose: bool = True) -> list[str]:
    """Greedy: try the FULL blacklist; if biomass grows, bisect to find the
    smallest subset that's actually needed for growth. Add those back."""
    # Try adding the full blacklist
    add_reactions_from(model, gsm, blacklist)
    bm = test_biomass(model, biomass_rxn, glc_uptake=glc_uptake)
    if verbose:
        print(f'biomass with FULL blacklist re-added: {bm:.6f}')
    if bm < biomass_threshold:
        # Not enough even with everything back — something else is missing
        remove_reactions(model, blacklist)
        if verbose:
            print('  biomass infeasible even with all blacklist; leaving them out')
        return []

    # Bisect: try without half, see if still grows
    keep = list(blacklist)
    remove_reactions(model, blacklist)

    needed = []
    pool = list(blacklist)
    while pool:
        candidate = pool[len(pool)//2:]   # try with right-half
        add_reactions_from(model, gsm, candidate)
        bm = test_biomass(model, biomass_rxn, glc_uptake=glc_uptake)
        if bm >= biomass_threshold:
            # Right half is enough — narrow further within that
            remove_reactions(model, candidate)
            if len(candidate) == 1:
                needed.append(candidate[0])
                add_reactions_from(model, gsm, candidate)
                pool = []
            else:
                pool = candidate
        else:
            # Need left half too
            remove_reactions(model, candidate)
            left = pool[:len(pool)//2]
            add_reactions_from(model, gsm, left + candidate)
            bm = test_biomass(model, biomass_rxn, glc_uptake=glc_uptake)
            if bm >= biomass_threshold:
                needed.extend(left)
                needed.extend(candidate)
            pool = []  # stop
    if verbose:
        print(f'  minimum blacklist subset needed for biomass: {len(needed)}')
    return needed
