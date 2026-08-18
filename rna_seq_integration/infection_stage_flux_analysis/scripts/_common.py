"""Shared constants and helpers for the infection-stage flux analysis (M1).

All of M1 reuses the *outputs* of the existing rna_seq_integration pipeline; no
FBA/MILP is re-run here. The single S1 transcriptome is a static expression prior
applied to every medium -- expression does NOT vary by infection stage, only flux
(via the medium) does. Keep that framing in mind for every metric here.
"""
from __future__ import annotations
import os
import re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# ---------------------------------------------------------------- paths
HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(HERE)                       # infection_stage_flux_analysis/
RSI = os.path.dirname(ANALYSIS)                        # rna_seq_integration/
OUTPUTS = os.path.join(RSI, "outputs")                 # frozen pipeline outputs (read-only)
SRC = os.path.join(RSI, "src")                         # importable pipeline modules

DATA = os.path.join(ANALYSIS, "data")
FLUX_RESULTS = os.path.join(ANALYSIS, "flux_results")
EXPR_RESULTS = os.path.join(ANALYSIS, "expression_results")
METHOD_CMP = os.path.join(ANALYSIS, "method_comparison")
PATHWAY = os.path.join(ANALYSIS, "pathway_analysis")
STATS = os.path.join(ANALYSIS, "statistics")

# External ModelSEED biochemistry DB (read-only source for KEGG map; snapshotted into data/)
MODELSEED_DB = ("/home/janakae/claude_projects/fungi_biochemistry_recon/"
                "data/external/ModelSEEDDatabase/Biochemistry")

# ---------------------------------------------------------------- constants
FLUX_EPS = 1e-6                       # matches stage1_overlay.py / stage2_*.py
METHODS = ["pfba", "eflux", "gimme", "imat"]
# Featured (headline) thresholds for the multi-threshold methods.
DEFAULT_THRESHOLD = {"gimme": "default", "imat": "default"}
# Major infection stages (cocktail = engineered mixes, reported separately).
MAJOR_STAGES = ["pre-infection", "biotrophic", "necrotrophic"]


def ensure_dirs():
    for d in (DATA, FLUX_RESULTS, EXPR_RESULTS, METHOD_CMP, PATHWAY, STATS):
        os.makedirs(d, exist_ok=True)


def strip_compartment(rxn_id: str) -> str:
    """rxn00248_m0 -> rxn00248 ; EX_cpd00027_e0 -> EX_cpd00027 (left as-is if no seed base)."""
    return re.sub(r"_[a-z][0-9]+$", "", rxn_id)


def seed_base(rxn_id: str):
    """Return the ModelSEED reaction base id (rxnNNNNN) or None for custom ids."""
    b = strip_compartment(rxn_id)
    return b if re.fullmatch(r"rxn[0-9]+", b) else None


def read_stage_map() -> pd.DataFrame:
    """condition_id, label, stage per condition -- taken from the frozen stage1 summary."""
    s = pd.read_csv(os.path.join(OUTPUTS, "stage1_summary.tsv"), sep="\t")
    m = s[["condition_id", "label", "stage"]].drop_duplicates("condition_id")
    return m.reset_index(drop=True)


def classify(flux, bin_, has_gpr):
    """Identical to stage1_overlay.classify -- agreement category per reaction."""
    nonzero = abs(flux) > FLUX_EPS
    if nonzero:
        if not has_gpr:
            return "ORPHAN_FLUX"
        if bin_ in ("hi", "med"):
            return "SUPPORTED"
        if bin_ == "lo":
            return "WEAK_SUPPORT"
        return "CONFLICT_FLUX_NO_EXPR"
    if bin_ == "hi":
        return "PRIMED_NOT_USED"
    return "SILENT_OK"


def _mcc(tp, tn, fp, fn):
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return (tp * tn - fp * fn) / denom if denom > 0 else 0.0


def concordance(df: pd.DataFrame, flux_col: str) -> dict:
    """Concordance of a method's flux against the S1 expression prior.

    Metrics (the first four replicate stage1_summary.tsv exactly when flux_col is
    the pFBA flux; the last two are new binary active-vs-expressed summaries):
      agreement_score, spearman_expr_vs_flux, hi_expr_recall, flux_precision_hi_med,
      mcc_active_vs_expr, jaccard_active_vs_expr
    Guards match stage1 (spearman -> 0 when <10 valid points or max|flux|~0).
    """
    per = df
    absf = per[flux_col].abs()
    has_gpr = per["n_genes"] > 0

    agree = [classify(f, b, hg) for f, b, hg
             in zip(per[flux_col], per["expression_bin"], has_gpr)]
    agree = pd.Series(agree, index=per.index)
    sup = int((agree == "SUPPORTED").sum())
    cfl = int((agree == "CONFLICT_FLUX_NO_EXPR").sum())
    pnu = int((agree == "PRIMED_NOT_USED").sum())
    denom = sup + cfl + pnu
    agreement_score = sup / denom if denom else 0.0

    gpr_mask = has_gpr
    expr_scores = per.loc[gpr_mask, "agg_mean_log2TPMp1"].to_numpy(dtype=float)
    abs_flux = absf[gpr_mask].to_numpy(dtype=float)
    valid = ~np.isnan(expr_scores)
    if valid.sum() > 10 and abs_flux[valid].max() > FLUX_EPS:
        rho, _ = spearmanr(expr_scores[valid], abs_flux[valid])
        rho = 0.0 if np.isnan(rho) else float(rho)
    else:
        rho = 0.0

    hi = per[per["expression_bin"] == "hi"]
    n_hi = len(hi)
    hi_recall = int((hi[flux_col].abs() > FLUX_EPS).sum()) / n_hi if n_hi else 0.0

    gpr_active = per[gpr_mask & (absf > FLUX_EPS)]
    n_active = len(gpr_active)
    flux_precision = (int(gpr_active["expression_bin"].isin(["hi", "med"]).sum())
                      / n_active) if n_active else 0.0

    # binary active (|flux|>eps) vs expressed (bin in hi/med), GPR'd reactions only
    g = per[gpr_mask]
    active = (g[flux_col].abs() > FLUX_EPS).to_numpy()
    expressed = g["expression_bin"].isin(["hi", "med"]).to_numpy()
    tp = int((active & expressed).sum()); tn = int((~active & ~expressed).sum())
    fp = int((active & ~expressed).sum()); fn = int((~active & expressed).sum())
    mcc = _mcc(tp, tn, fp, fn)
    jac = tp / (tp + fp + fn) if (tp + fp + fn) else 0.0

    return {
        "agreement_score": round(agreement_score, 4),
        "spearman_expr_vs_flux": round(rho, 4),
        "hi_expr_recall": round(hi_recall, 4),
        "flux_precision_hi_med": round(flux_precision, 4),
        "mcc_active_vs_expr": round(mcc, 4),
        "jaccard_active_vs_expr": round(jac, 4),
        "n_flux_nonzero": int((absf > FLUX_EPS).sum()),
        "n_gpr": int(gpr_mask.sum()),
    }


COMPOSITE_METRICS = ["spearman_expr_vs_flux", "hi_expr_recall",
                     "flux_precision_hi_med", "mcc_active_vs_expr"]


def add_composite(df: pd.DataFrame) -> pd.DataFrame:
    """Rank-normalize each component metric to [0,1] across all rows, mean = composite."""
    df = df.copy()
    for m in COMPOSITE_METRICS:
        r = df[m].rank(method="average")
        df[m + "_rn"] = (r - 1) / (len(df) - 1) if len(df) > 1 else 1.0
    df["composite"] = df[[m + "_rn" for m in COMPOSITE_METRICS]].mean(axis=1).round(4)
    return df
