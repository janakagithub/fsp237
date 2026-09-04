#!/opt/env/modelseed_cplex/bin/python3
"""Phase 4 — proteomics-driven context-specific flux: GIMME + iMAT.

Completes the flux-method panel for the proteome (Stage 1 gave E-Flux + the
enzyme-capacity model; this adds the two discrete context-extraction methods
that the RNA-seq track already has, so proteome and transcriptome are compared
on the same footing).

  * GIMME (Becker & Palsson 2008) — LP. bio >= 90% of max, minimise flux through
    reactions whose GPR-aggregated proteome abundance is below `lo_thr`.
  * iMAT  (Shlomi 2008) — MILP. bin reactions H/M/L on proteome abundance,
    maximise (#H carrying flux) + (#L inactive), bio >= 10% of max.

Expression = GPR-aggregated log2 TMT intensity per condition (max-of-OR /
min-of-AND), matching the RNA-seq convention. Thresholds are percentiles of the
*model-reaction* proteome-abundance distribution (recomputed for TMT; the
RNA-seq log2(TPM+1) cut-points do not transfer).

Reuses run_gimme / run_imat from the RNA-seq stage-2 module verbatim (same
solver, same MILP encoding), and the PDB media + condition map from phase1.

Outputs (proteomics_integration/outputs/context/):
  gimme_flux_matrix.tsv   rxn x condition GIMME flux (default threshold)
  imat_flux_matrix.tsv    rxn x condition iMAT flux (default threshold)
  imat_activity.tsv       per (condition, threshold) H_active / L_inactive tallies
  context_summary.tsv     one row per (condition, method, threshold)
  method_agreement.tsv    per reaction: active-call agreement across methods (PDA)
"""
import sys
from pathlib import Path

import cobra
import numpy as np
import pandas as pd

ROOT = Path("/home/janakae/fsp237")
SRC_RNA = ROOT / "rna_seq_integration/src"
SRC_PROT = ROOT / "proteomics_integration/src"
sys.path.insert(0, str(SRC_RNA))
sys.path.insert(0, str(SRC_PROT))
sys.path.insert(0, str(ROOT / "simulations"))

from gpr_expression import aggregate                       # noqa: E402
from run_simulation_panel import apply_media               # noqa: E402
from stage2_gimme_imat_eflux import run_gimme, run_imat    # noqa: E402
from phase1_eflux import PDB, CONDITIONS, BIO              # noqa: E402

MODEL_PATH = (ROOT / "simulations/gapfill_v1_v2/models/"
              "fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json")
MEANS = ROOT / "proteomics_integration/outputs/proteome_condition_means.tsv"
OUT = ROOT / "proteomics_integration/outputs/context"
OUT.mkdir(parents=True, exist_ok=True)
FLUX_EPS = 1e-6
ACT_EPS = 1e-3


def rxn_expr(model, gene_vals):
    """GPR-aggregated proteome abundance per reaction (log2 TMT)."""
    out = {}
    for r in model.reactions:
        s, _, _ = aggregate(r.gene_reaction_rule or "", gene_vals)
        out[r.id] = s
    return out


def main():
    model = cobra.io.load_json_model(str(MODEL_PATH))
    means = pd.read_csv(MEANS, sep="\t").set_index("gene")
    print(f"model {MODEL_PATH.name}: {len(model.reactions)} rxns; "
          f"proteome {len(means)} genes")

    # per-condition GPR-aggregated proteome abundance (log2 TMT)
    expr = {}
    for lab, col, _ in CONDITIONS:
        gv = means[col].dropna().to_dict()
        expr[lab] = rxn_expr(model, gv)

    # thresholds from the pooled model-reaction abundance distribution
    pool = np.array([v for lab, _, _ in CONDITIONS for v in expr[lab].values()
                     if v is not None and not np.isnan(v)])
    THR = {
        "default": (float(np.percentile(pool, 25)), float(np.percentile(pool, 75))),
        "strict":  (float(np.percentile(pool, 10)), float(np.percentile(pool, 90))),
    }
    print("thresholds (log2 TMT):",
          {k: (round(v[0], 2), round(v[1], 2)) for k, v in THR.items()})

    gimme_flux = pd.DataFrame(index=[r.id for r in model.reactions])
    gimme_flux["name"] = [r.name or "" for r in model.reactions]
    imat_flux = gimme_flux.copy()
    summary, act_rows = [], []

    for lab, col, media_key in CONDITIONS:
        eb = expr[lab]
        for thr_name, (lo, hi) in THR.items():
            # --- GIMME ---
            with model:
                apply_media(model, PDB[media_key], aerobic=True)
                model.objective = BIO
                fx, bio_max, pen = run_gimme(model, eb, lo)
            n_active = int(sum(1 for v in (fx.values if fx is not None else [])
                               if abs(v) >= ACT_EPS)) if fx is not None else 0
            summary.append(dict(condition=lab, method="GIMME", threshold=thr_name,
                                bio_max=round(bio_max or 0, 6),
                                penalty=round(pen, 3) if pen is not None else None,
                                n_active=n_active, n_H=None, n_L=None))
            if thr_name == "default" and fx is not None:
                gimme_flux[f"flux_{lab}"] = [float(fx.get(r.id, 0.0))
                                             for r in model.reactions]

            # --- iMAT ---
            with model:
                apply_media(model, PDB[media_key], aerobic=True)
                model.objective = BIO
                bmax = float(model.optimize().objective_value or 0.0)
                ifx, stats = run_imat(model, eb, lo, hi, bmax)
            summary.append(dict(condition=lab, method="iMAT", threshold=thr_name,
                                bio_max=round(bmax, 6), penalty=None,
                                n_active=int(sum(1 for v in (ifx.values if ifx is not None else [])
                                                 if abs(v) >= ACT_EPS)) if ifx is not None else 0,
                                n_H=stats.get("H_total"), n_L=stats.get("L_total")))
            act_rows.append(dict(condition=lab, threshold=thr_name, **stats))
            if thr_name == "default" and ifx is not None:
                imat_flux[f"flux_{lab}"] = [float(ifx.get(r.id, 0.0))
                                            for r in model.reactions]
            print(f"  {lab:9s} {thr_name:8s}  GIMME pen={pen}  "
                  f"iMAT H {stats.get('H_active')}/{stats.get('H_total')} "
                  f"L {stats.get('L_inactive')}/{stats.get('L_total')}")

    gimme_flux.to_csv(OUT / "gimme_flux_matrix.tsv", sep="\t")
    imat_flux.to_csv(OUT / "imat_flux_matrix.tsv", sep="\t")
    pd.DataFrame(act_rows).to_csv(OUT / "imat_activity.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(OUT / "context_summary.tsv", sep="\t", index=False)

    # --- cross-method active-call agreement on the PDA baseline ---
    # E-Flux flux from phase1, enzyme-capacity flux from phase2
    ef = pd.read_csv(ROOT / "proteomics_integration/outputs/eflux/eflux_flux_matrix.tsv",
                     sep="\t", index_col=0)
    ec = pd.read_csv(ROOT / "proteomics_integration/outputs/ec/ec_flux_matrix.tsv",
                     sep="\t", index_col=0)
    rows = []
    for r in model.reactions:
        calls = {
            "eflux":  abs(float(ef.get("eflux_PDA", pd.Series()).get(r.id, 0.0) or 0.0)) >= ACT_EPS,
            "gimme":  abs(float(gimme_flux.get("flux_PDA", pd.Series()).get(r.id, 0.0) or 0.0)) >= ACT_EPS,
            "imat":   abs(float(imat_flux.get("flux_PDA", pd.Series()).get(r.id, 0.0) or 0.0)) >= ACT_EPS,
            "ec":     abs(float(ec.get("flux_PDA", pd.Series()).get(r.id, 0.0) or 0.0)) >= ACT_EPS,
        }
        n_on = sum(calls.values())
        if n_on == 0:
            continue
        rows.append(dict(rxn_id=r.id, name=r.name or "", n_methods_active=n_on,
                         **{f"{k}_on": bool(v) for k, v in calls.items()}))
    pd.DataFrame(rows).sort_values("n_methods_active", ascending=False).to_csv(
        OUT / "method_agreement.tsv", sep="\t", index=False)

    print(f"\nwrote GIMME/iMAT matrices + summary to {OUT}")
    print(f"method_agreement rows: {len(rows)}")


if __name__ == "__main__":
    main()
