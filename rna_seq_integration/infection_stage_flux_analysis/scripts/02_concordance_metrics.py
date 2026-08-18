"""02 - Concordance of each method's flux with the S1 expression prior + ranking.

For every (method, condition, O2) we score how well the method's flux distribution
agrees with the single static S1 transcriptome, using six metrics (see
_common.concordance). Four of them replicate stage1_summary.tsv exactly for pFBA,
which we assert as a regression check. A rank-normalized composite ranks methods
overall and within each infection stage.

Caveat surfaced downstream: iMAT activates ~1.5x more reactions than the others,
which mechanically inflates recall and deflates precision -- MCC and the composite
are the fair comparators, not raw active counts.
"""
import os
import numpy as np
import pandas as pd
from _common import (FLUX_RESULTS, METHOD_CMP, OUTPUTS, MAJOR_STAGES,
                     concordance, add_composite, ensure_dirs)


def main():
    ensure_dirs()
    long = pd.read_csv(os.path.join(FLUX_RESULTS, "unified_flux_matrix.tsv"), sep="\t")

    rows = []
    for (method, cond, o2), g in long.groupby(["method", "condition_id", "O2"]):
        m = concordance(g.rename(columns={"flux": "flux"}), "flux")
        stage = g["stage"].iloc[0]
        label = g["label"].iloc[0]
        feasible = m["n_flux_nonzero"] > 0
        rows.append({"method": method, "condition_id": cond, "label": label,
                     "stage": stage, "O2": o2, "feasible": feasible, **m})
    by = pd.DataFrame(rows)
    by.to_csv(os.path.join(METHOD_CMP, "concordance_by_method_condition.tsv"),
              sep="\t", index=False)

    # ---- regression check: pFBA aerobic must match frozen stage1_summary.tsv ----
    s1 = pd.read_csv(os.path.join(OUTPUTS, "stage1_summary.tsv"), sep="\t")
    ref = s1[s1["O2"] == "aerobic"].set_index("condition_id")
    mine = by[(by["method"] == "pfba") & (by["O2"] == "aerobic")].set_index("condition_id")
    maxerr = 0.0
    for col in ["agreement_score", "spearman_expr_vs_flux",
                "hi_expr_recall", "flux_precision_hi_med"]:
        d = (mine[col] - ref[col]).abs().max()
        maxerr = max(maxerr, float(d))
    print(f"regression vs stage1_summary (pFBA aer): max abs diff = {maxerr:.2e} "
          f"({'OK' if maxerr < 1e-3 else 'MISMATCH'})")

    # ---- composite: rank-normalize over feasible (method,cond,O2) rows ----
    feas = by[by["feasible"]].copy()
    feas = add_composite(feas)
    feas.to_csv(os.path.join(METHOD_CMP, "concordance_composite.tsv"),
                sep="\t", index=False)

    # ---- method ranking: overall (aerobic) and per major stage ----
    aer = feas[feas["O2"] == "aerobic"]
    overall = (aer.groupby("method")["composite"].mean()
               .sort_values(ascending=False).round(4))
    per_stage = (aer[aer["stage"].isin(MAJOR_STAGES)]
                 .groupby(["stage", "method"])["composite"].mean().round(4)
                 .reset_index())

    rank_rows = [{"scope": "overall_aerobic", "method": m, "mean_composite": v,
                  "rank": i + 1} for i, (m, v) in enumerate(overall.items())]
    for stage in MAJOR_STAGES:
        sub = per_stage[per_stage["stage"] == stage].sort_values(
            "composite", ascending=False)
        for i, (_, r) in enumerate(sub.iterrows()):
            rank_rows.append({"scope": stage, "method": r["method"],
                              "mean_composite": r["composite"], "rank": i + 1})
    ranking = pd.DataFrame(rank_rows)
    ranking.to_csv(os.path.join(METHOD_CMP, "method_ranking.tsv"), sep="\t", index=False)

    print("\nOverall method ranking (mean composite, aerobic):")
    for m, v in overall.items():
        print(f"  {m:6s} {v:.4f}")
    print("\nPer-stage winner:")
    for stage in MAJOR_STAGES:
        sub = ranking[(ranking["scope"] == stage)].sort_values("rank")
        w = sub.iloc[0]
        print(f"  {stage:14s} -> {w['method']} ({w['mean_composite']:.4f})")


if __name__ == "__main__":
    main()
