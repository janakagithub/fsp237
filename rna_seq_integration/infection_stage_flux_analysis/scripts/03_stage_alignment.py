"""03 - Which infection-stage medium does the S1 transcriptome align with?

S1 is a single (pathogenic-state) transcriptome applied to every medium; only the
flux changes across media. So "alignment" asks: in which medium/stage does the
model's flux distribution best agree with the fixed S1 expression prior?

Two products:
  (A) stage_alignment.tsv        -- per (method, stage) mean composite concordance,
                                    ranked; plus the single best-aligned medium.
  (B) vanilla_fba_stage_view.tsv -- USER ADDITION: an explicit vanilla-FBA (pFBA)
                                    across-media comparison grouped by MAJOR stage,
                                    with the raw concordance components and the
                                    per-condition biomass, so the FBA-only picture
                                    is legible independent of the method comparison.
A stage is called a credible "best match" only if it tops the ranking under more
than one method (cross-method robustness table).
"""
import os
import pandas as pd
from _common import (METHOD_CMP, STATS, OUTPUTS, MAJOR_STAGES, ensure_dirs)


def main():
    ensure_dirs()
    comp = pd.read_csv(os.path.join(METHOD_CMP, "concordance_composite.tsv"), sep="\t")
    aer = comp[comp["O2"] == "aerobic"].copy()

    # (A) per (method, stage) alignment
    rows = []
    for (method, stage), g in aer.groupby(["method", "stage"]):
        best = g.sort_values("composite", ascending=False).iloc[0]
        rows.append({
            "method": method, "stage": stage,
            "n_conditions": len(g),
            "mean_composite": round(g["composite"].mean(), 4),
            "max_composite": round(g["composite"].max(), 4),
            "best_condition": best["condition_id"],
            "best_label": best["label"],
            "mean_spearman": round(g["spearman_expr_vs_flux"].mean(), 4),
            "mean_mcc": round(g["mcc_active_vs_expr"].mean(), 4),
        })
    align = pd.DataFrame(rows).sort_values(["method", "mean_composite"],
                                           ascending=[True, False])
    align.to_csv(os.path.join(STATS, "stage_alignment.tsv"), sep="\t", index=False)

    # cross-method robustness: rank of each MAJOR stage per method
    maj = align[align["stage"].isin(MAJOR_STAGES)].copy()
    maj["stage_rank"] = maj.groupby("method")["mean_composite"].rank(
        ascending=False, method="min").astype(int)
    robust = maj.pivot_table(index="stage", columns="method",
                             values="stage_rank", aggfunc="first")
    robust["times_ranked_1st"] = (robust == 1).sum(axis=1)
    robust = robust.sort_values("times_ranked_1st", ascending=False)
    robust.to_csv(os.path.join(STATS, "stage_alignment_robustness.tsv"), sep="\t")

    # single best-aligned medium overall (per method) + which stage it belongs to
    best_medium = (aer.sort_values("composite", ascending=False)
                   .groupby("method").head(1)
                   [["method", "condition_id", "label", "stage", "composite"]])
    best_medium.to_csv(os.path.join(STATS, "best_aligned_medium.tsv"),
                       sep="\t", index=False)

    # (B) vanilla-FBA (pFBA) across-media, grouped by MAJOR stage (user addition)
    s1 = pd.read_csv(os.path.join(OUTPUTS, "stage1_summary.tsv"), sep="\t")
    s1b = s1[s1["O2"] == "aerobic"][["condition_id", "biomass"]]
    fba = aer[aer["method"] == "pfba"].merge(s1b, on="condition_id", how="left")
    fba = fba[fba["stage"].isin(MAJOR_STAGES)].copy()
    fba_view = fba[["stage", "condition_id", "label", "biomass",
                    "spearman_expr_vs_flux", "hi_expr_recall",
                    "flux_precision_hi_med", "mcc_active_vs_expr",
                    "agreement_score", "composite"]].sort_values(
                        ["stage", "composite"], ascending=[True, False])
    fba_view.to_csv(os.path.join(STATS, "vanilla_fba_stage_view.tsv"),
                    sep="\t", index=False)
    fba_stage = (fba.groupby("stage")
                 .agg(n=("condition_id", "size"),
                      mean_biomass=("biomass", "mean"),
                      mean_composite=("composite", "mean"),
                      mean_spearman=("spearman_expr_vs_flux", "mean"))
                 .round(4).reindex(MAJOR_STAGES))
    fba_stage.to_csv(os.path.join(STATS, "vanilla_fba_stage_summary.tsv"), sep="\t")

    # ---- report ----
    print("Stage alignment (mean composite, aerobic) per method:")
    for method in ["pfba", "eflux", "gimme", "imat"]:
        sub = align[(align["method"] == method) &
                    (align["stage"].isin(MAJOR_STAGES))].sort_values(
                        "mean_composite", ascending=False)
        order = " > ".join(f"{r['stage']}({r['mean_composite']:.3f})"
                           for _, r in sub.iterrows())
        print(f"  {method:6s}: {order}")
    print("\nCross-method robustness (times a stage ranked #1 across methods):")
    print(robust["times_ranked_1st"].to_string())
    print("\nVanilla-FBA per-stage summary:")
    print(fba_stage.to_string())


if __name__ == "__main__":
    main()
