"""06 - Bundle M1 results into the website payload for the new top-level tab.

Writes a standalone atp-safe/infection_stage_payload.json consumed by an
ADDITIVE fetch in index.html. Nothing in the existing RNA-seq data flow
(reactions.json / rnaseq_payload.json) is touched.
"""
import os
import json
import pandas as pd
from _common import (METHOD_CMP, STATS, PATHWAY, DATA, ANALYSIS, MAJOR_STAGES)

SITE_JSON = "/home/janakae/fsp237/atp-safe/infection_stage_payload.json"


def rd(path, **kw):
    return pd.read_csv(path, sep="\t", **kw)


def main():
    ranking = rd(os.path.join(METHOD_CMP, "method_ranking.tsv"))
    comp = rd(os.path.join(METHOD_CMP, "concordance_composite.tsv"))
    align = rd(os.path.join(STATS, "stage_alignment.tsv"))
    robust = rd(os.path.join(STATS, "stage_alignment_robustness.tsv"))
    best = rd(os.path.join(STATS, "best_aligned_medium.tsv"))
    fba_view = rd(os.path.join(STATS, "vanilla_fba_stage_view.tsv"))
    fba_sum = rd(os.path.join(STATS, "vanilla_fba_stage_summary.tsv"))
    boot = rd(os.path.join(STATS, "stage_alignment_bootstrap.tsv"))
    with open(os.path.join(STATS, "bootstrap_summary.json")) as fh:
        boot_meta = json.load(fh)
    stage_diff = rd(os.path.join(PATHWAY, "stage_differential_pathways.tsv"))
    discord = rd(os.path.join(PATHWAY, "method_discordant_pathways.tsv"))
    kegg = rd(os.path.join(DATA, "rxn_kegg_map.tsv")).fillna("")

    aer = comp[comp["O2"] == "aerobic"].copy()

    # KEGG coverage
    n = len(kegg)
    kegg_cov = {
        "n_reactions": int(n),
        "kegg_reaction": int((kegg["kegg_reaction"].astype(str) != "").sum()),
        "kegg_pathway": int((kegg["kegg_pathway_id"].astype(str) != "").sum()),
        "distinct_kegg_pathways": len({p for c in kegg["kegg_pathway_id"]
                                       if isinstance(c, str) and c
                                       for p in c.split(";")}),
    }

    # winner (overall) + best-aligned stage (cross-method robustness)
    overall = ranking[ranking["scope"] == "overall_aerobic"].sort_values("rank")
    best_stage_row = robust.sort_values("times_ranked_1st", ascending=False).iloc[0]

    payload = {
        "meta": {
            "model": "FSP237 V10 (1622 rxns / 1268 mets / 1274 genes)",
            "dataset_label": "S1",
            "n_methods": 4,
            "n_conditions_aerobic": int(aer["condition_id"].nunique()),
            "featured_thresholds": {"gimme": "default", "imat": "default"},
            "best_method_overall": overall.iloc[0]["method"],
            "best_aligned_stage": str(best_stage_row["stage"]),
            "limitation": (
                "A single S1 transcriptome (technical replicates G1/G2/G3), "
                "representing the pathogenic state, is applied to every medium. "
                "Expression does NOT vary by infection stage — only flux (via "
                "the medium) does. So there is no differential-expression analysis; "
                "instead S1 is a fixed expression prior, and we ask which stage's "
                "media the model's flux best matches, and which integration method "
                "best reproduces S1."),
            "kegg_coverage": kegg_cov,
        },
        "method_ranking": ranking.to_dict("records"),
        "concordance": aer[[
            "method", "condition_id", "label", "stage",
            "spearman_expr_vs_flux", "hi_expr_recall", "flux_precision_hi_med",
            "mcc_active_vs_expr", "agreement_score", "composite",
            "n_flux_nonzero"]].round(4).to_dict("records"),
        "stage_alignment": align.to_dict("records"),
        "stage_robustness": robust.to_dict("records"),
        "best_medium": best.round(4).to_dict("records"),
        "vanilla_fba_view": fba_view.round(4).to_dict("records"),
        "vanilla_fba_summary": fba_sum.round(4).to_dict("records"),
        "bootstrap": {"B": boot_meta["B"], "note": boot_meta["note"],
                      "stages": boot.to_dict("records")},
        "stage_diff_pathways": stage_diff.round(4).head(15).to_dict("records"),
        "method_discordant_pathways": discord.round(4).head(15).to_dict("records"),
    }

    with open(SITE_JSON, "w") as fh:
        json.dump(payload, fh, indent=2)
    # also keep a copy inside the analysis folder's web/ for provenance
    web_dir = os.path.join(ANALYSIS, "web")
    os.makedirs(web_dir, exist_ok=True)
    with open(os.path.join(web_dir, "infection_stage_payload.json"), "w") as fh:
        json.dump(payload, fh, indent=2)

    print(f"wrote {SITE_JSON}")
    print(f"  best method overall : {payload['meta']['best_method_overall']}")
    print(f"  best-aligned stage  : {payload['meta']['best_aligned_stage']}")
    print(f"  KEGG coverage       : {kegg_cov}")


if __name__ == "__main__":
    main()
