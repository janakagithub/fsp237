#!/opt/env/modelseed_cplex/bin/python3
"""Stage 2 — GIMME + iMAT + E-Flux context-specific analyses.

Uses the CPLEX-backed environment because iMAT is a MILP with hundreds of
binary variables per condition. GLPK on the same problem takes 5+ minutes
per condition; CPLEX solves in seconds.


Three independent expression-integration techniques run per (condition, O2):

  * E-Flux (Colijn 2009) — continuous. Reaction bounds are scaled by
    (aggregated expression / 99th-%ile expression). Reactions without
    GPR keep default bounds. Then pFBA is solved on the constrained model.

  * GIMME (Becker & Palsson 2008) — LP. Enforce biomass ≥ 90 % of
    unconstrained max, minimize sum of |flux| in reactions whose
    aggregated expression falls below `lo_thr`. Reactions retained despite
    the penalty are the "must-carry" set.

  * iMAT (Shlomi 2008) — MILP. Reactions binned as
    H (expr ≥ hi_thr) / M / L (expr < lo_thr).
    Maximize (# H reactions that carry flux) + (# L reactions inactive).

Threshold pairs on Mean log2(TPM+1) over model genes:
  default : (lo=25 %ile, hi=75 %ile)   = (3.12, 6.46)
  strict  : (lo=10 %ile, hi=90 %ile)   = (0.88, 7.97)
  narrow  : (lo=50 %ile, hi=75 %ile)   = (4.99, 6.46)

Outputs:
  outputs/stage2_eflux/<cond>_<o2>.tsv
  outputs/stage2_gimme/<cond>_<o2>_<thr>.tsv
  outputs/stage2_imat/<cond>_<o2>_<thr>.tsv
  outputs/stage2_summary.tsv                # one row per (cond, O2, method, threshold)
"""
import json
import sys
from pathlib import Path

import cobra
import numpy as np
import pandas as pd

sys.path.insert(0, '/home/janakae/fsp237/simulations')
from run_simulation_panel import CONDITIONS, apply_media  # noqa: E402

ROOT = Path('/home/janakae/fungalTemplate/imm904CobraModel/rna_seq_integration')
MODEL_PATH = Path('/home/janakae/fsp237/simulations/gapfill_v1_v2/models/'
                  'fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json')
REXP_TSV = ROOT / 'outputs' / 'reaction_expression.tsv'
COV_JSON = ROOT / 'outputs' / 'coverage_summary.json'
SUMMARY_OUT = ROOT / 'outputs' / 'stage2_summary.tsv'

FLUX_EPS = 1e-6
GIMME_BIO_FRAC = 0.90
IMAT_BIO_FRAC = 0.10          # biomass floor for iMAT (biological viability)
IMAT_ACTIVITY_EPS = 1e-3      # flux magnitude that counts a reaction "on"
IMAT_MAX_FLUX = 1000.0

THRESHOLD_PAIRS = {
    'default':  (3.121, 6.457),   # 25%–75% of model-gene Mean log2(TPM+1)
    'strict':   (0.882, 7.966),   # 10%–90%
    # 'narrow' (lo=median) was tried but produced an iMAT MILP too large for
    # per-batch runtime (~343 L reactions × binary indicators). Kept only for
    # GIMME (LP, tractable).
    'narrow':   (4.993, 6.457),   # 50%–75%
}
IMAT_THRESHOLD_PAIRS = {k: v for k, v in THRESHOLD_PAIRS.items() if k != 'narrow'}
IMAT_TIME_LIMIT_SEC = 60     # per-solve wallclock cap; accept best incumbent
IMAT_MIP_GAP = 0.05          # 5% optimality gap is fine for our reporting


# ---------------------------------------------------------------------------

def load_expression():
    rexp = pd.read_csv(REXP_TSV, sep='\t', comment='#').set_index('rxn_id')
    return rexp


# ---------------------------------------------------------------------------
# E-Flux
# ---------------------------------------------------------------------------

def run_eflux(model, expr_by_rxn):
    """Scale reaction bounds by expression / P99(expression)."""
    valid = [v for v in expr_by_rxn.values()
             if v is not None and not (isinstance(v, float) and np.isnan(v))]
    if not valid:
        return None
    p99 = float(np.quantile(valid, 0.99))
    if p99 <= 0:
        return None

    orig_bounds = {}
    for r in model.reactions:
        expr = expr_by_rxn.get(r.id)
        if expr is None or (isinstance(expr, float) and np.isnan(expr)):
            continue
        scale = min(1.0, expr / p99)
        orig_bounds[r.id] = (r.lower_bound, r.upper_bound)
        if r.upper_bound > 0:
            r.upper_bound = r.upper_bound * scale
        if r.lower_bound < 0:
            r.lower_bound = r.lower_bound * scale
    try:
        sol = cobra.flux_analysis.pfba(model)
    except Exception:
        sol = model.optimize()
    # Restore
    for rid, (lb, ub) in orig_bounds.items():
        r = model.reactions.get_by_id(rid)
        r.lower_bound, r.upper_bound = lb, ub
    return sol


# ---------------------------------------------------------------------------
# GIMME
# ---------------------------------------------------------------------------

def run_gimme(model, expr_by_rxn, lo_thr):
    """LP: enforce bio ≥ 90 % max, minimize sum of |flux| * (lo_thr - expr) for
    reactions with expr < lo_thr. Returns (fluxes, bio_max, penalty)."""
    bio = model.reactions.get_by_id('bio_gsm')
    sol_max = model.optimize()
    if sol_max.status != 'optimal' or sol_max.objective_value < FLUX_EPS:
        return None, sol_max.objective_value, None
    bio_max = sol_max.objective_value
    bio.lower_bound = GIMME_BIO_FRAC * bio_max

    problem = model.problem
    cons_vars = []
    obj_terms = []
    for r in model.reactions:
        expr = expr_by_rxn.get(r.id)
        if expr is None or (isinstance(expr, float) and np.isnan(expr)):
            continue
        if expr >= lo_thr:
            continue
        penalty = float(lo_thr - expr)
        v_pos = problem.Variable(f'_gm_p_{r.id}', lb=0, ub=IMAT_MAX_FLUX)
        v_neg = problem.Variable(f'_gm_n_{r.id}', lb=0, ub=IMAT_MAX_FLUX)
        c = problem.Constraint(r.flux_expression - v_pos + v_neg, lb=0, ub=0,
                                name=f'_gm_abs_{r.id}')
        cons_vars.extend([v_pos, v_neg, c])
        obj_terms.append(penalty * (v_pos + v_neg))
    if not cons_vars:
        return sol_max.fluxes, bio_max, 0.0

    model.add_cons_vars(cons_vars)
    model.objective = problem.Objective(sum(obj_terms), direction='min')
    sol = model.optimize()
    fluxes = sol.fluxes if sol.status == 'optimal' else None
    penalty = sol.objective_value if sol.status == 'optimal' else None
    model.remove_cons_vars(cons_vars)
    return fluxes, bio_max, penalty


# ---------------------------------------------------------------------------
# iMAT
# ---------------------------------------------------------------------------

def run_imat(model, expr_by_rxn, lo_thr, hi_thr, bio_max):
    """MILP: max (#H active) + (#L inactive), keep bio ≥ IMAT_BIO_FRAC * bio_max."""
    if bio_max < FLUX_EPS:
        return None, {'H_active': 0, 'H_total': 0, 'L_inactive': 0, 'L_total': 0}
    bio = model.reactions.get_by_id('bio_gsm')
    bio.lower_bound = IMAT_BIO_FRAC * bio_max

    problem = model.problem
    cons_vars = []
    obj_terms = []
    h_total = l_total = 0

    for r in model.reactions:
        expr = expr_by_rxn.get(r.id)
        if expr is None or (isinstance(expr, float) and np.isnan(expr)):
            continue
        if expr >= hi_thr:
            h_total += 1
            # active if |v| >= eps; two indicator binaries y_pos, y_neg summing into y
            # y = 1 iff v >= eps or v <= -eps. big-M linking.
            y_pos = problem.Variable(f'_im_yp_{r.id}', type='binary')
            y_neg = problem.Variable(f'_im_yn_{r.id}', type='binary')
            # v >= eps - M*(1 - y_pos)   ->  v - eps + M*(1-y_pos) >= 0
            #                                 v + M - M*y_pos >= eps
            c1 = problem.Constraint(r.flux_expression + IMAT_MAX_FLUX * (1 - y_pos),
                                     lb=IMAT_ACTIVITY_EPS, name=f'_im_hp_{r.id}')
            # v <= -eps + M*(1 - y_neg)
            c2 = problem.Constraint(r.flux_expression - IMAT_MAX_FLUX * (1 - y_neg),
                                     ub=-IMAT_ACTIVITY_EPS, name=f'_im_hn_{r.id}')
            cons_vars.extend([y_pos, y_neg, c1, c2])
            obj_terms.append(y_pos + y_neg)
        elif expr < lo_thr:
            l_total += 1
            # inactive if -eps <= v <= eps. z=1 iff inactive.
            z = problem.Variable(f'_im_z_{r.id}', type='binary')
            # v <=  eps + M*(1 - z)
            c1 = problem.Constraint(r.flux_expression - IMAT_MAX_FLUX * (1 - z),
                                     ub=IMAT_ACTIVITY_EPS, name=f'_im_lu_{r.id}')
            # v >= -eps - M*(1 - z)
            c2 = problem.Constraint(r.flux_expression + IMAT_MAX_FLUX * (1 - z),
                                     lb=-IMAT_ACTIVITY_EPS, name=f'_im_ll_{r.id}')
            cons_vars.extend([z, c1, c2])
            obj_terms.append(z)

    if not obj_terms:
        return None, {'H_active': 0, 'H_total': 0, 'L_inactive': 0, 'L_total': 0}

    model.add_cons_vars(cons_vars)
    model.objective = problem.Objective(sum(obj_terms), direction='max')
    # CPLEX-specific tuning: time + gap limits keep the MILP tractable.
    try:
        model.solver.configuration.timeout = IMAT_TIME_LIMIT_SEC
        prob = model.solver.problem
        prob.parameters.mip.tolerances.mipgap.set(IMAT_MIP_GAP)
        prob.parameters.timelimit.set(IMAT_TIME_LIMIT_SEC)
    except Exception:
        pass
    sol = model.optimize()
    if sol.status not in ('optimal', 'feasible', 'time_limit'):
        model.remove_cons_vars(cons_vars)
        return None, {'H_active': 0, 'H_total': h_total, 'L_inactive': 0, 'L_total': l_total}

    fluxes = sol.fluxes
    h_act = sum(1 for r in model.reactions
                 if (expr_by_rxn.get(r.id) is not None
                     and not (isinstance(expr_by_rxn.get(r.id), float) and np.isnan(expr_by_rxn.get(r.id)))
                     and expr_by_rxn[r.id] >= hi_thr
                     and abs(fluxes.get(r.id, 0)) >= IMAT_ACTIVITY_EPS))
    l_inact = sum(1 for r in model.reactions
                   if (expr_by_rxn.get(r.id) is not None
                       and not (isinstance(expr_by_rxn.get(r.id), float) and np.isnan(expr_by_rxn.get(r.id)))
                       and expr_by_rxn[r.id] < lo_thr
                       and abs(fluxes.get(r.id, 0)) < IMAT_ACTIVITY_EPS))
    model.remove_cons_vars(cons_vars)
    return fluxes, {'H_active': h_act, 'H_total': h_total,
                    'L_inactive': l_inact, 'L_total': l_total,
                    'objective': int(sol.objective_value)}


# ---------------------------------------------------------------------------

def dump_flux_table(rexp, fluxes, out_path, extra_cols=None):
    df = rexp.copy()
    df['flux'] = df.index.map(lambda rid: float(fluxes.get(rid, 0.0)))
    df['is_active'] = df['flux'].abs() > FLUX_EPS
    if extra_cols:
        for k, v in extra_cols.items():
            df[k] = v
    df.reset_index().to_csv(out_path, sep='\t', index=False)


def main():
    print(f'model  : {MODEL_PATH.name}')
    print(f'rexp   : {REXP_TSV.name}')
    model = cobra.io.load_json_model(str(MODEL_PATH))
    rexp = load_expression()
    expr_by_rxn = rexp['agg_mean_log2TPMp1'].to_dict()

    for d in ('stage2_eflux', 'stage2_gimme', 'stage2_imat'):
        (ROOT / 'outputs' / d).mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for cond_id, label, stage, c_sources, notes in CONDITIONS:
        for aerobic in (True, False):
            o2_tag = 'aerobic' if aerobic else 'anaerobic'
            print(f'\n=== {cond_id}  {o2_tag} ===')

            # --- baseline pFBA to get bio_max (used by GIMME + iMAT) ---
            with model:
                apply_media(model, c_sources, aerobic)
                for r in model.reactions:
                    r.objective_coefficient = 1 if r.id == 'bio_gsm' else 0
                sol0 = model.optimize()
                bio_max = sol0.objective_value if sol0.status == 'optimal' else 0.0
                if bio_max is None or bio_max < FLUX_EPS:
                    print(f'  infeasible baseline (bio={bio_max}) — skipping method runs')
                    for method in ('eflux',):
                        summary_rows.append({
                            'condition_id': cond_id, 'label': label, 'stage': stage,
                            'O2': o2_tag, 'method': method, 'threshold': 'n/a',
                            'status': 'infeasible', 'biomass': 0.0, 'n_active': 0,
                            'n_H_active': 0, 'n_H_total': 0,
                            'n_L_inactive': 0, 'n_L_total': 0, 'objective': None,
                        })
                    continue

                # --- E-FLUX (no threshold; just continuous scaling) ---
                with model:
                    sol_e = run_eflux(model, expr_by_rxn)
                if sol_e is not None and sol_e.status == 'optimal':
                    bio_e = float(sol_e.fluxes.get('bio_gsm', 0.0))
                    n_act = int((sol_e.fluxes.abs() > FLUX_EPS).sum())
                    dump_flux_table(rexp, sol_e.fluxes,
                                     ROOT / 'outputs' / 'stage2_eflux' / f'{cond_id}_{o2_tag}.tsv',
                                     extra_cols={'method': 'eflux'})
                    print(f'  eflux                     : bio={bio_e:.4f}   active={n_act}')
                    summary_rows.append({
                        'condition_id': cond_id, 'label': label, 'stage': stage,
                        'O2': o2_tag, 'method': 'eflux', 'threshold': 'p99',
                        'status': 'optimal', 'biomass': round(bio_e, 6),
                        'n_active': n_act, 'n_H_active': None, 'n_H_total': None,
                        'n_L_inactive': None, 'n_L_total': None, 'objective': None,
                    })
                else:
                    print(f'  eflux                     : infeasible')
                    summary_rows.append({
                        'condition_id': cond_id, 'label': label, 'stage': stage,
                        'O2': o2_tag, 'method': 'eflux', 'threshold': 'p99',
                        'status': 'infeasible', 'biomass': 0.0, 'n_active': 0,
                        'n_H_active': None, 'n_H_total': None,
                        'n_L_inactive': None, 'n_L_total': None, 'objective': None,
                    })

                # --- GIMME (threshold sweep) ---
                for thr_name, (lo, hi) in THRESHOLD_PAIRS.items():
                    with model:
                        fluxes, bio_max_g, penalty = run_gimme(model, expr_by_rxn, lo)
                    if fluxes is not None:
                        n_act = int((fluxes.abs() > FLUX_EPS).sum())
                        bio_g = float(fluxes.get('bio_gsm', 0.0))
                        dump_flux_table(rexp, fluxes,
                                         ROOT / 'outputs' / 'stage2_gimme' /
                                         f'{cond_id}_{o2_tag}_{thr_name}.tsv',
                                         extra_cols={'method': 'gimme',
                                                     'threshold': thr_name})
                        print(f'  gimme lo={lo:5.2f} ({thr_name:<7s}): '
                              f'bio={bio_g:.4f}  active={n_act}  penalty={penalty:.3f}')
                        summary_rows.append({
                            'condition_id': cond_id, 'label': label, 'stage': stage,
                            'O2': o2_tag, 'method': 'gimme', 'threshold': thr_name,
                            'status': 'optimal', 'biomass': round(bio_g, 6),
                            'n_active': n_act, 'n_H_active': None, 'n_H_total': None,
                            'n_L_inactive': None, 'n_L_total': None,
                            'objective': round(penalty, 3) if penalty is not None else None,
                        })
                    else:
                        summary_rows.append({
                            'condition_id': cond_id, 'label': label, 'stage': stage,
                            'O2': o2_tag, 'method': 'gimme', 'threshold': thr_name,
                            'status': 'infeasible', 'biomass': 0.0, 'n_active': 0,
                            'n_H_active': None, 'n_H_total': None,
                            'n_L_inactive': None, 'n_L_total': None, 'objective': None,
                        })

                # --- iMAT (threshold sweep, MILP time-limited) ---
                for thr_name, (lo, hi) in IMAT_THRESHOLD_PAIRS.items():
                    with model:
                        fluxes, stats = run_imat(model, expr_by_rxn, lo, hi, bio_max)
                    if fluxes is not None:
                        bio_i = float(fluxes.get('bio_gsm', 0.0))
                        n_act = int((fluxes.abs() > FLUX_EPS).sum())
                        dump_flux_table(rexp, fluxes,
                                         ROOT / 'outputs' / 'stage2_imat' /
                                         f'{cond_id}_{o2_tag}_{thr_name}.tsv',
                                         extra_cols={'method': 'imat',
                                                     'threshold': thr_name})
                        print(f'  imat  ({thr_name:<7s})         : '
                              f'bio={bio_i:.4f}  active={n_act}  '
                              f'H={stats["H_active"]}/{stats["H_total"]}  '
                              f'L={stats["L_inactive"]}/{stats["L_total"]}  obj={stats.get("objective")}')
                        summary_rows.append({
                            'condition_id': cond_id, 'label': label, 'stage': stage,
                            'O2': o2_tag, 'method': 'imat', 'threshold': thr_name,
                            'status': 'optimal', 'biomass': round(bio_i, 6),
                            'n_active': n_act,
                            'n_H_active': stats['H_active'], 'n_H_total': stats['H_total'],
                            'n_L_inactive': stats['L_inactive'], 'n_L_total': stats['L_total'],
                            'objective': stats.get('objective'),
                        })
                    else:
                        print(f'  imat  ({thr_name:<7s})         : infeasible')
                        summary_rows.append({
                            'condition_id': cond_id, 'label': label, 'stage': stage,
                            'O2': o2_tag, 'method': 'imat', 'threshold': thr_name,
                            'status': 'infeasible', 'biomass': 0.0, 'n_active': 0,
                            'n_H_active': stats['H_active'], 'n_H_total': stats['H_total'],
                            'n_L_inactive': stats['L_inactive'], 'n_L_total': stats['L_total'],
                            'objective': None,
                        })

    sdf = pd.DataFrame(summary_rows)
    sdf.to_csv(SUMMARY_OUT, sep='\t', index=False)
    print(f'\nwrote {SUMMARY_OUT}')

    # Compact final summary
    print('\n=== compact roll-up (aerobic, default threshold) ===')
    aer = sdf[(sdf.O2 == 'aerobic') & (sdf.threshold.isin(['default', 'p99']))].copy()
    print(aer[['condition_id', 'stage', 'method', 'biomass', 'n_active',
                'n_H_active', 'n_H_total', 'n_L_inactive', 'n_L_total']]
          .to_string(index=False))


if __name__ == '__main__':
    main()
