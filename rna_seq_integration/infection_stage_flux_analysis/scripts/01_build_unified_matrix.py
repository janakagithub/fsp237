"""01 - Build the unified reaction x method x condition x O2 flux matrix.

Reshapes the frozen per-condition TSVs from the four methods (vanilla pFBA,
E-Flux, GIMME, iMAT) into one tidy long table, merged with the S1 expression
prior (method-invariant), the 23-bucket pathway assignment, the KEGG map, and
the condition->stage tag. GIMME/iMAT are filtered to the featured DEFAULT
threshold for the headline table; a companion all-thresholds table is also written.

No FBA is re-run; agreement is recomputed per method with the frozen classify().
Infeasible conditions (empty files / all-zero flux, e.g. several anaerobic FA and
pentose media) are tolerated: their rows carry flux 0 and are flagged downstream.
"""
import os
import glob
import re
import pandas as pd
from _common import (OUTPUTS, DATA, FLUX_RESULTS, DEFAULT_THRESHOLD, FLUX_EPS,
                     classify, read_stage_map, ensure_dirs)

# leading columns shared by every per-condition file
KEEP_META = ["rxn_id", "name", "gpr", "n_genes", "n_genes_with_expr",
             "agg_mean_log2TPMp1", "agg_mean_TPM", "expression_bin"]

FNAME_RE = re.compile(r"^(?P<cond>.+?)_(?P<o2>aerobic|anaerobic)(?:_(?P<thr>\w+))?\.tsv$")


def _parse_name(path):
    m = FNAME_RE.match(os.path.basename(path))
    return m.group("cond"), m.group("o2"), m.group("thr")


def _load_method(subdir, flux_col, method, only_threshold=None):
    """Return a long df for one method across all its per-condition files."""
    rows = []
    for path in sorted(glob.glob(os.path.join(OUTPUTS, subdir, "*.tsv"))):
        cond, o2, thr = _parse_name(path)
        if only_threshold is not None and thr != only_threshold:
            continue
        try:
            df = pd.read_csv(path, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        if df.empty or flux_col not in df.columns:
            continue
        d = df[[c for c in KEEP_META if c in df.columns]].copy()
        d["flux"] = df[flux_col].astype(float)
        d["method"] = method
        d["condition_id"] = cond
        d["O2"] = o2
        d["threshold"] = thr if thr else "default"
        rows.append(d)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def build(all_thresholds: bool) -> pd.DataFrame:
    specs = [
        ("stage1_overlay", "flux_pFBA", "pfba", None),
        ("stage2_eflux", "flux", "eflux", None),
        ("stage2_gimme", "flux", "gimme",
         None if all_thresholds else DEFAULT_THRESHOLD["gimme"]),
        ("stage2_imat", "flux", "imat",
         None if all_thresholds else DEFAULT_THRESHOLD["imat"]),
    ]
    parts = [_load_method(sd, fc, m, thr) for sd, fc, m, thr in specs]
    long = pd.concat([p for p in parts if not p.empty], ignore_index=True)

    long["abs_flux"] = long["flux"].abs()
    long["is_active"] = long["abs_flux"] > FLUX_EPS
    long["has_gpr"] = long["n_genes"] > 0
    long["agreement"] = [classify(f, b, hg) for f, b, hg
                         in zip(long["flux"], long["expression_bin"], long["has_gpr"])]

    # static per-reaction attributes
    kegg = pd.read_csv(os.path.join(DATA, "rxn_kegg_map.tsv"), sep="\t").fillna("")
    pa = pd.read_csv(os.path.join(OUTPUTS, "pathway_assignment.tsv"), sep="\t")
    kegg = kegg.merge(pa[["rxn_id", "compartment_bucket"]], on="rxn_id", how="left")
    long = long.merge(
        kegg[["rxn_id", "pathway_bucket", "compartment_bucket", "kegg_reaction",
              "kegg_pathway_id", "kegg_pathway_name", "ec_number"]],
        on="rxn_id", how="left").rename(columns={"pathway_bucket": "pathway"})

    # condition -> stage tag
    long = long.merge(read_stage_map()[["condition_id", "stage", "label"]],
                      on="condition_id", how="left")
    return long


def main():
    ensure_dirs()

    featured = build(all_thresholds=False)
    cols = ["rxn_id", "method", "condition_id", "O2", "stage", "label",
            "flux", "abs_flux", "is_active", "name", "gpr", "n_genes",
            "n_genes_with_expr", "agg_mean_log2TPMp1", "agg_mean_TPM",
            "expression_bin", "has_gpr", "agreement", "pathway",
            "compartment_bucket", "kegg_reaction", "kegg_pathway_id",
            "kegg_pathway_name", "ec_number"]
    featured = featured[cols]
    fp = os.path.join(FLUX_RESULTS, "unified_flux_matrix.tsv")
    featured.to_csv(fp, sep="\t", index=False)

    all_thr = build(all_thresholds=True)
    all_thr = all_thr[["threshold"] + [c for c in cols if c in all_thr.columns]]
    ap = os.path.join(FLUX_RESULTS, "unified_flux_matrix_allthresholds.tsv")
    all_thr.to_csv(ap, sep="\t", index=False)

    # wide flux pivots (rxn x condition) per method, aerobic -- for Escher/heatmaps
    aer = featured[featured["O2"] == "aerobic"]
    for method in featured["method"].unique():
        w = (aer[aer["method"] == method]
             .pivot_table(index="rxn_id", columns="condition_id",
                          values="flux", aggfunc="first"))
        w.to_csv(os.path.join(FLUX_RESULTS, f"wide_flux_{method}_aerobic.tsv"), sep="\t")

    # report
    print(f"wrote {fp}")
    print(f"  rows ............. {len(featured)}")
    combos = featured.groupby("method")[["condition_id", "O2"]].apply(
        lambda d: d.drop_duplicates().shape[0])
    for m, n in combos.items():
        print(f"  {m:6s} (cond,O2) combos = {n}")
    print(f"  distinct reactions {featured['rxn_id'].nunique()}")
    print(f"wrote {ap}  rows {len(all_thr)}")

    # integrity: expression columns must be identical across methods per reaction
    chk = (featured.groupby("rxn_id")["agg_mean_log2TPMp1"].nunique() > 1).sum()
    print(f"  integrity: reactions with method-varying expression = {chk} (expect 0)")


if __name__ == "__main__":
    main()
