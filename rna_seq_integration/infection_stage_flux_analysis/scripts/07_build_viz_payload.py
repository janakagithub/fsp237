"""07 - Build the M2 visualization payload (separate, lazy-loaded).

M1 shipped a compact infection_stage_payload.json for the Overview / Method
Comparison / Stage Alignment / Pathways sub-tabs. M2 adds heavier visual
sub-tabs (pathway heatmap, infection-stage dashboard, flux-expression
agreement, reaction/gene table). To keep the M1 tab fast, that heavier data
goes into a SEPARATE atp-safe/infection_stage_viz.json that index.html fetches
only when an M2 sub-tab is first opened.

Nothing here re-runs FBA/MILP; it reshapes the frozen unified matrix and the
pathway matrices already on disk. Aerobic media only (the featured panel).
"""
import os
import json
import numpy as np
import pandas as pd
from _common import (FLUX_RESULTS, PATHWAY, ANALYSIS, METHODS, MAJOR_STAGES,
                     FLUX_EPS, ensure_dirs)

SITE_JSON = "/home/janakae/fsp237/atp-safe/infection_stage_viz.json"
AGREEMENT_ORDER = ["SUPPORTED", "WEAK_SUPPORT", "PRIMED_NOT_USED",
                   "SILENT_OK", "ORPHAN_FLUX", "CONFLICT_FLUX_NO_EXPR"]
DISCORDANT = {"PRIMED_NOT_USED", "ORPHAN_FLUX", "CONFLICT_FLUX_NO_EXPR"}
KEGG_MIN_RXN = 3   # KEGG pathways with >=3 model reactions enter the heatmap


def heatmap_block(long, group_col, min_rxn=1):
    """Nested {pathways:[{name,n_rxn}], cells:{method:{cond:{pathway:rate}}}}."""
    sub = long[long[group_col].notna() & (long[group_col].astype(str) != "")].copy()
    # n_rxn is constant per pathway (reaction set doesn't depend on condition)
    nrxn = (sub.groupby(group_col)["rxn_id"].nunique()
            .sort_values(ascending=False))
    keep = nrxn[nrxn >= min_rxn].index.tolist()
    sub = sub[sub[group_col].isin(keep)]
    g = (sub.groupby([group_col, "method", "condition_id"])
         .agg(n=("rxn_id", "nunique"), a=("is_active", "sum")).reset_index())
    g["rate"] = (g["a"] / g["n"]).round(3)
    cells = {m: {} for m in METHODS}
    for (pw, m, c), r in zip(zip(g[group_col], g["method"], g["condition_id"]),
                             g["rate"]):
        cells.setdefault(m, {}).setdefault(c, {})[pw] = float(r)
    pathways = [{"name": p, "n_rxn": int(nrxn[p])} for p in keep]
    return {"pathways": pathways, "cells": cells,
            "n_omitted": int((nrxn < min_rxn).sum())}


def main():
    ensure_dirs()
    long = pd.read_csv(os.path.join(FLUX_RESULTS, "unified_flux_matrix.tsv"),
                       sep="\t")
    aer = long[long["O2"] == "aerobic"].copy()

    # conditions (id -> label, stage), ordered by id
    cond = (aer[["condition_id", "label", "stage"]].drop_duplicates()
            .sort_values("condition_id"))
    conditions = cond.to_dict("records")
    cond_ids = cond["condition_id"].tolist()

    # media grouped by major/other stage (for the dashboard)
    stage_media = {}
    for r in conditions:
        stage_media.setdefault(r["stage"], []).append(
            {"condition_id": r["condition_id"], "label": r["label"]})

    # ---- heatmaps (two taxonomies) ----
    heat = {
        "bucket": heatmap_block(aer, "pathway", min_rxn=1),
        "kegg": heatmap_block(aer, "kegg_pathway_name", min_rxn=KEGG_MIN_RXN),
    }

    # ---- agreement summary (per method x condition, GPR'd reactions) ----
    gpr = aer[aer["has_gpr"] == True]  # noqa: E712
    ag = (gpr.groupby(["method", "condition_id", "label", "stage", "agreement"])
          .size().reset_index(name="n"))
    summary = []
    for (m, c, lab, st), grp in ag.groupby(["method", "condition_id",
                                            "label", "stage"]):
        cats = {k: 0 for k in AGREEMENT_ORDER}
        for _, row in grp.iterrows():
            cats[row["agreement"]] = int(row["n"])
        summary.append({"method": m, "condition_id": c, "label": lab,
                        "stage": st, "cats": cats, "total": int(grp["n"].sum())})

    # ---- reaction table (one row per reaction; per-method aggregates) ----
    # static (expression-side) fields are method-invariant -> take from pfba slice
    base = (aer[aer["method"] == "pfba"]
            .groupby("rxn_id")
            .agg(name=("name", "first"), pathway=("pathway", "first"),
                 kegg_reaction=("kegg_reaction", "first"),
                 kegg_pathway_name=("kegg_pathway_name", "first"),
                 ec_number=("ec_number", "first"),
                 expression_bin=("expression_bin", "first"),
                 has_gpr=("has_gpr", "first"),
                 log2tpm=("agg_mean_log2TPMp1", "first"),
                 n_genes=("n_genes", "first")).reset_index())

    # per (rxn, method): active fraction, mean |flux|, dominant agreement
    def dominant(s):
        return s.value_counts().idxmax() if len(s) else ""
    perm = (aer.groupby(["rxn_id", "method"])
            .agg(af=("is_active", "mean"),
                 mf=("abs_flux", "mean"),
                 ag=("agreement", dominant),
                 disc=("agreement", lambda s: float(np.mean([x in DISCORDANT
                                                             for x in s])))).reset_index())

    perm_idx = {(r.rxn_id, r.method): r for r in perm.itertuples(index=False)}
    rows = []
    for b in base.itertuples(index=False):
        rec = {
            "rxn": b.rxn_id, "name": b.name or "",
            "path": b.pathway or "", "kegg": b.kegg_reaction or "",
            "kpath": b.kegg_pathway_name or "", "ec": b.ec_number or "",
            "bin": b.expression_bin or "absent", "gpr": bool(b.has_gpr),
            "log2tpm": (round(float(b.log2tpm), 3)
                        if pd.notna(b.log2tpm) else None),
            "ngenes": int(b.n_genes) if pd.notna(b.n_genes) else 0,
            "m": {},
        }
        for m in METHODS:
            pr = perm_idx.get((b.rxn_id, m))
            if pr is not None:
                rec["m"][m] = {"af": round(float(pr.af), 3),
                               "mf": round(float(pr.mf), 4),
                               "ag": pr.ag,
                               "disc": round(float(pr.disc), 3)}
            else:
                rec["m"][m] = {"af": 0.0, "mf": 0.0, "ag": "", "disc": 0.0}
        rows.append(rec)

    payload = {
        "meta": {
            "note": ("M2 visualization payload (lazy-loaded). Aerobic media only; "
                     "GIMME/iMAT at the default threshold. Heatmap KEGG taxonomy "
                     "shows pathways with >= %d model reactions." % KEGG_MIN_RXN),
            "methods": METHODS,
            "conditions": conditions,
            "condition_order": cond_ids,
            "major_stages": MAJOR_STAGES,
            "agreement_order": AGREEMENT_ORDER,
            "agreement_glossary": {
                "SUPPORTED": "active flux & high/med expression",
                "WEAK_SUPPORT": "active flux & low expression",
                "PRIMED_NOT_USED": "expressed (hi/med) but carries no flux",
                "SILENT_OK": "no expression & no flux (correctly off)",
                "ORPHAN_FLUX": "carries flux but no/absent expression",
                "CONFLICT_FLUX_NO_EXPR": "flux with an explicit no-expression call",
            },
        },
        "stage_media": stage_media,
        "heatmap": heat,
        "agreement": {"order": AGREEMENT_ORDER, "summary": summary},
        "reaction_table": rows,
    }

    with open(SITE_JSON, "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    web_dir = os.path.join(ANALYSIS, "web")
    os.makedirs(web_dir, exist_ok=True)
    with open(os.path.join(web_dir, "infection_stage_viz.json"), "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))

    sz = os.path.getsize(SITE_JSON) / 1024
    print(f"wrote {SITE_JSON} ({sz:.0f} KB)")
    print(f"  heatmap bucket pathways : {len(heat['bucket']['pathways'])}")
    print(f"  heatmap KEGG pathways   : {len(heat['kegg']['pathways'])} "
          f"(omitted {heat['kegg']['n_omitted']} with <{KEGG_MIN_RXN} rxns)")
    print(f"  agreement rows          : {len(summary)}")
    print(f"  reaction_table rows     : {len(rows)}")


if __name__ == "__main__":
    main()
