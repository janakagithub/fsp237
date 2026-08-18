"""05 - Pathway x method x condition activity, in both taxonomies.

Aggregates the unified matrix to the pathway level so we can see (a) which
pathways are most stage-differential in FLUX and (b) where the four methods most
DISAGREE about whether a pathway is active. Two taxonomies are emitted:
  - the primary 23-bucket classifier ('pathway')
  - the KEGG pathway grouping ('kegg_pathway_name'; ~41% of reactions carry one)

For each (pathway, method, condition, O2) we report active_rate = fraction of the
pathway's reactions carrying |flux|>eps, and total_flux = sum|flux|. Downstream
rankings use aerobic rows.
"""
import os
import numpy as np
import pandas as pd
from _common import (FLUX_RESULTS, PATHWAY, MAJOR_STAGES, FLUX_EPS, ensure_dirs)


def _agg(long, group_col):
    sub = long[long[group_col].notna() & (long[group_col].astype(str) != "")]
    g = (sub.groupby([group_col, "method", "condition_id", "O2", "stage"])
         .agg(n_rxn=("rxn_id", "nunique"),
              n_active=("is_active", "sum"),
              total_flux=("abs_flux", "sum"))
         .reset_index())
    g["active_rate"] = (g["n_active"] / g["n_rxn"]).round(4)
    g["total_flux"] = g["total_flux"].round(4)
    return g


def main():
    ensure_dirs()
    long = pd.read_csv(os.path.join(FLUX_RESULTS, "unified_flux_matrix.tsv"), sep="\t")

    bucket = _agg(long, "pathway")
    bucket.to_csv(os.path.join(PATHWAY, "pathway_method_matrix.tsv"),
                  sep="\t", index=False)
    kegg = _agg(long, "kegg_pathway_name")
    kegg.to_csv(os.path.join(PATHWAY, "kegg_pathway_matrix.tsv"),
                sep="\t", index=False)

    aer = bucket[bucket["O2"] == "aerobic"]

    # (a) most stage-differential pathways (by pFBA active_rate spread across MAJOR stages)
    pf = aer[(aer["method"] == "pfba") & (aer["stage"].isin(MAJOR_STAGES))]
    stage_mean = (pf.groupby(["pathway", "stage"])["active_rate"].mean()
                  .unstack("stage").reindex(columns=MAJOR_STAGES))
    stage_mean["range"] = (stage_mean.max(axis=1) - stage_mean.min(axis=1)).round(4)
    stage_diff = stage_mean.sort_values("range", ascending=False).round(4)
    stage_diff.to_csv(os.path.join(PATHWAY, "stage_differential_pathways.tsv"), sep="\t")

    # (b) most method-discordant pathways (spread of active_rate across methods, per condition, averaged)
    piv = (aer.pivot_table(index=["pathway", "condition_id"], columns="method",
                           values="active_rate", aggfunc="first"))
    piv["method_range"] = piv.max(axis=1) - piv.min(axis=1)
    discord = (piv.groupby("pathway")["method_range"].mean()
               .sort_values(ascending=False).round(4))
    discord.to_frame("mean_method_range").to_csv(
        os.path.join(PATHWAY, "method_discordant_pathways.tsv"), sep="\t")

    print(f"wrote pathway_method_matrix.tsv ({len(bucket)} rows), "
          f"kegg_pathway_matrix.tsv ({len(kegg)} rows)")
    print("\nTop 8 stage-differential pathways (pFBA active-rate range across major stages):")
    print(stage_diff.head(8).to_string())
    print("\nTop 8 method-discordant pathways (mean active-rate range across methods):")
    print(discord.head(8).to_string())


if __name__ == "__main__":
    main()
