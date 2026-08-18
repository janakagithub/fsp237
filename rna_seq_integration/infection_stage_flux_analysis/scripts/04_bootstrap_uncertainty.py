"""04 - Gene-level bootstrap uncertainty on the S1 -> stage alignment.

Only ONE biological transcriptome exists (S1; G1/G2/G3 are technical replicates),
so no p-value on stage differences is valid. Instead we quantify how sensitive the
alignment is to the particular set of measured genes via a NON-PARAMETRIC bootstrap
over model genes (NOT over the 3 technical replicates): resample model genes with
replacement, re-aggregate every GPR, re-bin on the resampled gene distribution, and
recompute the stage ranking.

The bootstrap uses the VANILLA pFBA flux, which is independent of expression, so
resampling expression is self-consistent. (E-Flux/GIMME/iMAT fluxes were computed
FROM the point-estimate expression; bootstrapping them would require re-solving
1000x and is out of scope for M1.)

Outputs (descriptive only): statistics/stage_alignment_bootstrap.tsv and
statistics/bootstrap_summary.json (per-stage composite intervals + fraction of
replicates in which each stage ranks #1).
"""
import os
import sys
import json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from _common import (OUTPUTS, FLUX_RESULTS, STATS, MAJOR_STAGES, FLUX_EPS,
                     concordance, COMPOSITE_METRICS, ensure_dirs)

sys.path.insert(0, os.path.join(os.path.dirname(OUTPUTS), "src"))
from gpr_expression import parse_gpr  # exact GPR parser used by the pipeline

XLSX = "/home/janakae/fsp237/expression-data/S1_normalized_expression.xlsx"
B = 500
SEED = 1234


def fast_aggregate(clauses, expr):
    """score = max_OR(min_AND(expr[g])); nan if no listed gene has expression."""
    if not clauses:
        return np.nan
    or_maxes = []
    for cl in clauses:
        vals = [expr[g] for g in cl if g in expr]
        if vals:
            or_maxes.append(min(vals))
    return max(or_maxes) if or_maxes else np.nan


def bin_scores(agg, cuts_lo, cuts_hi):
    out = np.empty(len(agg), dtype=object)
    for i, s in enumerate(agg):
        if s is None or (isinstance(s, float) and np.isnan(s)) or s <= 0:
            out[i] = "absent"
        elif s >= cuts_hi:
            out[i] = "hi"
        elif s >= cuts_lo:
            out[i] = "med"
        else:
            out[i] = "lo"
    return out


def main():
    ensure_dirs()
    rng = np.random.default_rng(SEED)

    # reactions + parsed GPRs
    rexp = pd.read_csv(os.path.join(OUTPUTS, "reaction_expression.tsv"),
                       sep="\t", comment="#")
    rexp = rexp[["rxn_id", "gpr", "n_genes"]].copy()
    clauses_by_rxn = [parse_gpr(g if isinstance(g, str) else "") for g in rexp["gpr"]]

    # gene -> mean log2(TPM+1)
    xl = pd.read_excel(XLSX, sheet_name="Normalized log2 TPM")
    mean_by_gene = dict(zip(xl["Gene ID"], xl["Mean log2(TPM+1)"].astype(float)))
    model_genes = sorted({g for cls in clauses_by_rxn for cl in cls for g in cl})
    covered = [g for g in model_genes
               if g in mean_by_gene and not np.isnan(mean_by_gene[g])]
    covered = np.array(covered)

    # fixed vanilla pFBA aerobic flux per condition
    long = pd.read_csv(os.path.join(FLUX_RESULTS, "unified_flux_matrix.tsv"), sep="\t")
    pf = long[(long["method"] == "pfba") & (long["O2"] == "aerobic")]
    stage_of = dict(zip(pf["condition_id"], pf["stage"]))
    conditions = sorted(pf["condition_id"].unique())
    flux_by_cond = {c: dict(zip(pf[pf["condition_id"] == c]["rxn_id"],
                               pf[pf["condition_id"] == c]["flux"]))
                    for c in conditions}
    n_genes_arr = rexp["n_genes"].to_numpy()

    def one_replicate(expr):
        # re-aggregate + re-bin
        agg = np.array([fast_aggregate(cls, expr) for cls in clauses_by_rxn], dtype=float)
        gene_scores = np.array([expr[g] for g in expr], dtype=float)
        lo, hi = np.percentile(gene_scores, [25, 75])
        bins = bin_scores(agg, lo, hi)
        base = pd.DataFrame({"rxn_id": rexp["rxn_id"], "n_genes": n_genes_arr,
                             "agg_mean_log2TPMp1": agg, "expression_bin": bins})
        # concordance per condition
        recs = []
        for c in conditions:
            df = base.copy()
            df["flux"] = df["rxn_id"].map(flux_by_cond[c]).fillna(0.0)
            m = concordance(df, "flux")
            recs.append({"condition_id": c, "stage": stage_of[c], **m})
        cc = pd.DataFrame(recs)
        # composite = rank-normalized mean across conditions (within replicate)
        for m in COMPOSITE_METRICS:
            r = cc[m].rank(method="average")
            cc[m + "_rn"] = (r - 1) / (len(cc) - 1) if len(cc) > 1 else 1.0
        cc["composite"] = cc[[m + "_rn" for m in COMPOSITE_METRICS]].mean(axis=1)
        maj = cc[cc["stage"].isin(MAJOR_STAGES)]
        return maj.groupby("stage")["composite"].mean()

    # point estimate
    expr0 = {g: mean_by_gene[g] for g in covered}
    point = one_replicate(expr0).reindex(MAJOR_STAGES)

    # bootstrap
    boot = {s: [] for s in MAJOR_STAGES}
    first_counts = {s: 0 for s in MAJOR_STAGES}
    for b in range(B):
        drawn = rng.choice(covered, size=len(covered), replace=True)
        # duplicates don't change min/max; the resampled *distribution* shifts bins
        expr = {g: mean_by_gene[g] for g in np.unique(drawn)}
        sm = one_replicate(expr).reindex(MAJOR_STAGES)
        for s in MAJOR_STAGES:
            boot[s].append(float(sm[s]))
        first_counts[sm.idxmax()] += 1

    rows, summ = [], {"B": B, "seed": SEED,
                      "note": ("gene-level bootstrap on expression-independent pFBA "
                               "flux; descriptive intervals only, no p-values "
                               "(single biological transcriptome)")}
    for s in MAJOR_STAGES:
        arr = np.array(boot[s])
        lo, hi = np.percentile(arr, [2.5, 97.5])
        rows.append({"stage": s,
                     "point_composite": round(float(point[s]), 4),
                     "boot_mean": round(float(arr.mean()), 4),
                     "ci95_lo": round(float(lo), 4),
                     "ci95_hi": round(float(hi), 4),
                     "frac_ranked_1st": round(first_counts[s] / B, 4)})
    out = pd.DataFrame(rows).sort_values("boot_mean", ascending=False)
    out.to_csv(os.path.join(STATS, "stage_alignment_bootstrap.tsv"),
               sep="\t", index=False)
    summ["stages"] = rows
    with open(os.path.join(STATS, "bootstrap_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=2)

    print(f"gene bootstrap B={B} on pFBA (expression-independent flux)")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
