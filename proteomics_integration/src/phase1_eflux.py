#!/opt/env/modelseed_cplex/bin/python3
"""Phase 1 — proteomics-driven E-Flux on condition-matched PDB media.

The MPLEx proteome (PDA / half / onetenth) is matched 1:1 to the simulated PDB
dilution media (19/20/21). Unlike the single S1 transcriptome (one state applied
to every medium), here EXPRESSION VARIES BY CONDITION, so E-Flux gives a genuine
per-condition flux contrast on top of the media dilution.

E-Flux (Colijn 2009): reaction bounds are scaled by (enzyme abundance / P99),
capped at 1. We use LINEAR intensity (2**log2) — the faithful Colijn formulation
and the right handling for relative TMT (positive, proportional to protein amount;
the compressed log2 range would barely constrain anything). GPR aggregation is the
standard fungal convention max_OR(min_AND) reused from gpr_expression.py.

For every condition we also solve an UNCONSTRAINED pFBA baseline on the same
medium, so each reaction's E-Flux effect is measured against its own medium
(isolating the proteome contribution from the dilution contribution).

Outputs (proteomics_integration/outputs/eflux/):
  <cond>_eflux.tsv          per-reaction: baseline vs eflux flux, scale, expr
  eflux_flux_matrix.tsv     rxn x {baseline,eflux} flux across the 3 conditions
  eflux_summary.tsv         one row per condition (biomass, #active, #bound)
  eflux_de_concordance.tsv  reaction flux-change vs proteome DE logFC
  eflux_report.md           narrative + numbers for the website tab
"""
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd
import cobra

ROOT = Path("/home/janakae/fsp237")
sys.path.insert(0, str(ROOT / "rna_seq_integration" / "src"))
sys.path.insert(0, str(ROOT / "simulations"))
from gpr_expression import aggregate                    # noqa: E402
from run_simulation_panel import apply_media            # noqa: E402

MODEL = ROOT / "simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json"
MEANS = ROOT / "proteomics_integration/outputs/proteome_condition_means.tsv"
OUT   = ROOT / "proteomics_integration/outputs/eflux"
OUT.mkdir(parents=True, exist_ok=True)
BIO = "bio_gsm"
FLUX_EPS = 1e-6

# condition -> (mean column, matched PDB media c_sources)
from importlib.util import spec_from_file_location, module_from_spec
# PDB c_sources live in the conditions snippet; inline them to stay self-contained.
PDB = {
 "19_pdb_baseline": {'cpd00027_e0':5,'cpd00082_e0':0.00917529,'cpd00076_e0':0.0265603,'cpd00035_e0':0.00565912,'cpd00051_e0':0.0125258,'cpd00132_e0':0.0679374,'cpd00041_e0':0.0158347,'cpd00084_e0':0.0033426,'cpd00023_e0':0.0163471,'cpd00053_e0':0.0461466,'cpd00033_e0':0.00539483,'cpd00119_e0':0.00324936,'cpd00322_e0':0.00485144,'cpd00107_e0':0.00321329,'cpd00039_e0':0.0050883,'cpd00060_e0':0.00382207,'cpd00066_e0':0.00455308,'cpd00129_e0':0.00610209,'cpd00054_e0':0.00613453,'cpd00161_e0':0.00457939,'cpd00065_e0':0.00198301,'cpd00069_e0':0.00278255,'cpd00156_e0':0.0186256,'cpd00305_e0':2.72483e-05,'cpd00220_e0':7.3884e-06,'cpd00218_e0':0.000744052,'cpd00644_e0':0.000118939,'cpd00215_e0':0.000153423,'cpd00393_e0':3.0242e-06},
 "20_pdb_half": {'cpd00027_e0':2.5,'cpd00082_e0':0.00458764,'cpd00076_e0':0.0132802,'cpd00035_e0':0.00282956,'cpd00051_e0':0.00626288,'cpd00132_e0':0.0339687,'cpd00041_e0':0.00791735,'cpd00084_e0':0.0016713,'cpd00023_e0':0.00817353,'cpd00053_e0':0.0230733,'cpd00033_e0':0.00269741,'cpd00119_e0':0.00162468,'cpd00322_e0':0.00242572,'cpd00107_e0':0.00160665,'cpd00039_e0':0.00254415,'cpd00060_e0':0.00191104,'cpd00066_e0':0.00227654,'cpd00129_e0':0.00305104,'cpd00054_e0':0.00306726,'cpd00161_e0':0.00228969,'cpd00065_e0':0.000991504,'cpd00069_e0':0.00139128,'cpd00156_e0':0.00931279,'cpd00305_e0':1.36242e-05,'cpd00220_e0':3.6942e-06,'cpd00218_e0':0.000372026,'cpd00644_e0':5.94695e-05,'cpd00215_e0':7.67113e-05,'cpd00393_e0':1.5121e-06},
 "21_pdb_onetenth": {'cpd00027_e0':0.5,'cpd00082_e0':0.000917529,'cpd00076_e0':0.00265603,'cpd00035_e0':0.000565912,'cpd00051_e0':0.00125257,'cpd00132_e0':0.00679375,'cpd00041_e0':0.00158347,'cpd00084_e0':0.00033426,'cpd00023_e0':0.00163471,'cpd00053_e0':0.00461466,'cpd00033_e0':0.000539483,'cpd00119_e0':0.000324936,'cpd00322_e0':0.000485145,'cpd00107_e0':0.000321329,'cpd00039_e0':0.00050883,'cpd00060_e0':0.000382207,'cpd00066_e0':0.000455308,'cpd00129_e0':0.000610209,'cpd00054_e0':0.000613453,'cpd00161_e0':0.000457939,'cpd00065_e0':0.000198301,'cpd00069_e0':0.000278255,'cpd00156_e0':0.00186256,'cpd00305_e0':2.7248e-06,'cpd00220_e0':7.388e-07,'cpd00218_e0':7.44052e-05,'cpd00644_e0':1.18939e-05,'cpd00215_e0':1.53423e-05,'cpd00393_e0':3.024e-07},
}
CONDITIONS = [   # (label, mean_col, media_key)
    ("PDA",      "PDA_mean",      "19_pdb_baseline"),
    ("half",     "half_mean",     "20_pdb_half"),
    ("onetenth", "onetenth_mean", "21_pdb_onetenth"),
]


def eflux_scale(model, expr_by_rxn):
    """Return dict rxn_id -> scale in [0,1] and the p99 used."""
    vals = [v for v in expr_by_rxn.values() if v is not None and np.isfinite(v)]
    p99 = float(np.quantile(vals, 0.99)) if vals else 0.0
    scale = {}
    if p99 > 0:
        for rid, v in expr_by_rxn.items():
            if v is not None and np.isfinite(v):
                scale[rid] = min(1.0, v / p99)
    return scale, p99


def apply_eflux(model, scale):
    orig = {}
    for r in model.reactions:
        s = scale.get(r.id)
        if s is None:
            continue
        orig[r.id] = (r.lower_bound, r.upper_bound)
        if r.upper_bound > 0:
            r.upper_bound *= s
        if r.lower_bound < 0:
            r.lower_bound *= s
    return orig


def main():
    model = cobra.io.load_json_model(str(MODEL))
    print(f"model: {MODEL.name}  ({len(model.reactions)} rxns, {len(model.genes)} genes)")
    means = pd.read_csv(MEANS, sep="\t").set_index("gene")

    rxn_meta = {r.id: (r.name or "", r.subsystem or "", r.gene_reaction_rule or "")
                for r in model.reactions}
    flux_matrix = pd.DataFrame(index=[r.id for r in model.reactions])
    flux_matrix["name"] = [rxn_meta[r][0] for r in flux_matrix.index]
    flux_matrix["subsystem"] = [rxn_meta[r][1] for r in flux_matrix.index]

    summary = []
    per_cond_expr = {}
    for label, col, media_key in CONDITIONS:
        # linear intensity per gene
        lin = {g: float(2.0 ** v) for g, v in means[col].dropna().items()}
        expr_by_rxn = {}
        for r in model.reactions:
            score, ntot, nexpr = aggregate(r.gene_reaction_rule or "", lin)
            expr_by_rxn[r.id] = score if (nexpr > 0 and np.isfinite(score)) else np.nan
        per_cond_expr[label] = expr_by_rxn
        scale, p99 = eflux_scale(model, expr_by_rxn)
        n_scaled = len(scale)
        n_bound = sum(1 for s in scale.values() if s < 0.5)  # meaningfully throttled

        with model:
            apply_media(model, PDB[media_key], aerobic=True)
            model.objective = BIO
            # --- baseline pFBA (no expression constraint) ---
            base = cobra.flux_analysis.pfba(model)
            base_bio = float(base.fluxes.get(BIO, 0.0))
            base_active = int((base.fluxes.abs() > FLUX_EPS).sum())
            # --- E-Flux on same medium ---
            with model:
                apply_eflux(model, scale)
                try:
                    ef = cobra.flux_analysis.pfba(model)
                    ef_status = ef.status
                except Exception:
                    ef = model.optimize(); ef_status = ef.status
                ef_bio = float(ef.fluxes.get(BIO, 0.0))
                ef_active = int((ef.fluxes.abs() > FLUX_EPS).sum())

        # per-reaction table
        df = pd.DataFrame(index=[r.id for r in model.reactions])
        df["name"] = df.index.map(lambda i: rxn_meta[i][0])
        df["subsystem"] = df.index.map(lambda i: rxn_meta[i][1])
        df["gpr"] = df.index.map(lambda i: rxn_meta[i][2])
        df["expr_linear"] = df.index.map(lambda i: expr_by_rxn.get(i, np.nan))
        df["eflux_scale"] = df.index.map(lambda i: scale.get(i, np.nan))
        df["flux_baseline"] = df.index.map(lambda i: float(base.fluxes.get(i, 0.0)))
        df["flux_eflux"] = df.index.map(lambda i: float(ef.fluxes.get(i, 0.0)))
        df["dflux"] = df["flux_eflux"] - df["flux_baseline"]
        df.reset_index(names="rxn_id").to_csv(OUT / f"{label}_eflux.tsv", sep="\t", index=False)

        flux_matrix[f"base_{label}"] = df["flux_baseline"]
        flux_matrix[f"eflux_{label}"] = df["flux_eflux"]

        summary.append(dict(condition=label, media_key=media_key,
                            glucose_uptake=PDB[media_key]["cpd00027_e0"],
                            baseline_biomass=round(base_bio, 6),
                            eflux_biomass=round(ef_bio, 6),
                            eflux_status=ef_status,
                            biomass_ratio=round(ef_bio / base_bio, 4) if base_bio else None,
                            n_active_baseline=base_active, n_active_eflux=ef_active,
                            n_rxn_scaled=n_scaled, n_rxn_throttled=n_bound,
                            p99_linear=p99))
        print(f"  {label:8s} media={media_key:16s} base_bio={base_bio:.4f} "
              f"eflux_bio={ef_bio:.4f} ratio={ef_bio/base_bio if base_bio else 0:.3f} "
              f"active {base_active}->{ef_active} throttled={n_bound}/{n_scaled}")

    flux_matrix.to_csv(OUT / "eflux_flux_matrix.tsv", sep="\t")
    sdf = pd.DataFrame(summary)
    sdf.to_csv(OUT / "eflux_summary.tsv", sep="\t", index=False)

    # ---------------- capacity titration ----------------
    # Standard E-Flux is non-binding here (biomass headroom). Titrate a global
    # capacity multiplier C: reaction UB = C * expr_norm (expr_norm = linear/P99,
    # cap 1). As C shrinks, enzyme capacity tightens; the C at which biomass
    # first falls below its own-medium baseline is the proteome's capacity
    # headroom, and the reaction pinned at its ceiling is the first bottleneck.
    # Because the proteome differs by condition, the curves and bottlenecks
    # differ — this is the genuinely proteomics-driven readout.
    CAPS = [1000, 300, 100, 30, 10, 3, 1, 0.3, 0.1]
    tit_rows = []
    bottleneck_rows = []
    for label, col, media_key in CONDITIONS:
        expr_by_rxn = per_cond_expr[label]
        scale, p99 = eflux_scale(model, expr_by_rxn)   # expr_norm in [0,1]
        base_bio = float(sdf.loc[sdf.condition == label, "baseline_biomass"].iloc[0])
        for C in CAPS:
            with model:
                apply_media(model, PDB[media_key], aerobic=True)
                model.objective = BIO
                for r in model.reactions:
                    s = scale.get(r.id)
                    if s is None:
                        continue
                    ceil = C * s
                    if r.upper_bound > 0:
                        r.upper_bound = min(r.upper_bound, ceil)
                    if r.lower_bound < 0:
                        r.lower_bound = max(r.lower_bound, -ceil)
                sol = model.optimize()
                bio = float(sol.objective_value or 0.0) if sol.status == "optimal" else 0.0
                tit_rows.append(dict(condition=label, capacity=C,
                                     biomass=round(bio, 6),
                                     frac_of_baseline=round(bio / base_bio, 4) if base_bio else None))
                # bottleneck reactions at this C: GPR reaction pinned at its ceiling
                if sol.status == "optimal" and bio > FLUX_EPS and bio < 0.999 * base_bio:
                    for r in model.reactions:
                        s = scale.get(r.id)
                        if s is None or not (r.gene_reaction_rule or "").strip():
                            continue
                        f = float(sol.fluxes.get(r.id, 0.0))
                        ceil = C * s
                        if ceil > FLUX_EPS and abs(f) >= 0.98 * ceil:
                            bottleneck_rows.append(dict(condition=label, capacity=C,
                                rxn_id=r.id, name=r.name or "",
                                subsystem=r.subsystem or "", flux=round(f, 5),
                                ceiling=round(ceil, 5), expr_norm=round(s, 4)))
    pd.DataFrame(tit_rows).to_csv(OUT / "eflux_titration.tsv", sep="\t", index=False)
    bdf = pd.DataFrame(bottleneck_rows)
    bdf.to_csv(OUT / "eflux_bottlenecks.tsv", sep="\t", index=False)
    # first (largest C) bottleneck per condition
    first_bottleneck = {}
    if not bdf.empty:
        for label in bdf.condition.unique():
            sub = bdf[bdf.condition == label]
            Cmax = sub.capacity.max()
            first_bottleneck[label] = sub[sub.capacity == Cmax][
                ["rxn_id", "name", "subsystem", "capacity"]].to_dict("records")
    print("titration done; first bottlenecks:",
          {k: [r["rxn_id"] for r in v] for k, v in first_bottleneck.items()})

    # ---------------- DE concordance ----------------
    # For each contrast, does the E-Flux flux move in the same direction as the
    # proteome logFC of the reaction's controlling enzyme?
    conc_rows = []
    for contrast, (a, b) in {"half_vs_PDA": ("half", "PDA"),
                             "onetenth_vs_PDA": ("onetenth", "PDA"),
                             "onetenth_vs_half": ("onetenth", "half")}.items():
        fc_col = f"logFC_{contrast}"
        adj_col = f"adjP_{contrast}"
        if fc_col not in means.columns:
            continue
        for r in model.reactions:
            gpr = r.gene_reaction_rule or ""
            if not gpr.strip():
                continue
            # reaction-level proteome logFC = max_OR(min_AND) of gene logFC
            fc, ntot, nexpr = aggregate(gpr, means[fc_col].dropna().to_dict())
            if not np.isfinite(fc):
                continue
            fa = flux_matrix.at[r.id, f"eflux_{a}"]
            fb = flux_matrix.at[r.id, f"eflux_{b}"]
            if abs(fa) < FLUX_EPS and abs(fb) < FLUX_EPS:
                continue
            dflux = fa - fb
            conc_rows.append(dict(contrast=contrast, rxn_id=r.id,
                                  name=r.name or "", subsystem=r.subsystem or "",
                                  proteome_logFC=fc, flux_a=fa, flux_b=fb,
                                  dflux=dflux,
                                  concordant=int(np.sign(fc) == np.sign(dflux) and dflux != 0)))
    conc = pd.DataFrame(conc_rows)
    conc.to_csv(OUT / "eflux_de_concordance.tsv", sep="\t", index=False)

    # concordance rate per contrast (reactions with |dflux|>eps and |logFC|>0.5)
    conc_stats = {}
    if not conc.empty:
        for contrast in conc.contrast.unique():
            sub = conc[(conc.contrast == contrast) & (conc.dflux.abs() > FLUX_EPS)
                       & (conc.proteome_logFC.abs() > 0.5)]
            if len(sub):
                conc_stats[contrast] = dict(n=len(sub),
                                            concordant=int(sub.concordant.sum()),
                                            rate=round(sub.concordant.mean(), 3))

    # ---------------- report ----------------
    rep = [f"# Phase 1 — proteomics-driven E-Flux\n",
           f"Model {MODEL.name}; E-Flux (Colijn 2009) with LINEAR TMT intensity, "
           f"per-condition P99 normalization; GPR max_OR(min_AND). Baseline = "
           f"unconstrained pFBA on the same medium.\n",
           "## Per-condition summary\n",
           "```\n" + sdf.to_string(index=False) + "\n```\n",
           "## Capacity titration (proteome-implied bottleneck onset)\n",
           "Biomass as a fraction of own-medium baseline as the global capacity "
           "multiplier C shrinks. First bottleneck reaction(s) per condition: "
           + json.dumps({k: [r["rxn_id"] for r in v] for k, v in first_bottleneck.items()})
           + "\n",
           "## E-Flux ↔ proteome DE concordance\n",
           "Reactions with |Δflux|>0 and |proteome logFC|>0.5; fraction whose "
           "E-Flux flux change agrees in sign with the proteome fold-change.\n"]
    for c, s in conc_stats.items():
        rep.append(f"- **{c}**: {s['concordant']}/{s['n']} = {s['rate']*100:.1f}% concordant")
    rep.append("\n## Outputs\n- `<cond>_eflux.tsv`, `eflux_flux_matrix.tsv`, "
               "`eflux_summary.tsv`, `eflux_de_concordance.tsv`")
    (OUT / "eflux_report.md").write_text("\n".join(rep))
    (OUT / "eflux_concordance_stats.json").write_text(json.dumps(conc_stats, indent=2))
    print("\nconcordance:", conc_stats)
    print("wrote outputs to", OUT)


if __name__ == "__main__":
    main()
