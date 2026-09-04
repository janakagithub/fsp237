#!/opt/env/modelseed_cplex/bin/python3
"""Phase 2 — proteome-allocated enzyme-capacity model (sMOMENT/GECKO-style).

E-Flux (Phase 1) caps flux by raw abundance; that ignores that a slow, heavy
enzyme buys less flux per gram than a fast, light one. Here we convert measured
abundance into *capacity* the proper way:

    enzyme cost   c_i = MW_i / kcat_i           (g protein per unit flux)
    capacity      v_i <= share_i(cond) * B * kcat_i / MW_i
                      =  a_i(cond) / c_i

where share_i(cond) is the reaction's fraction of the measured proteome mass in
that condition (GPR-aggregated relative TMT, linearized), and B is a single
enzyme-budget knob swept over a wide range. Because the proteome differs by
condition, the per-reaction ceilings — and hence the flux state and the
saturated (budget-limiting) enzymes — differ by condition. This is the genuine
proteomics-driven flux prediction.

kcat is an EC-class prior (Bar-Even 2011 median ~10 s^-1, refined by EC first
digit); MW is computed from the FSP237 proteome sequences. Both are
ORDER-OF-MAGNITUDE and swappable for BRENDA/SABIO/DLKcat kcats later without
touching the framework.

Outputs (proteomics_integration/outputs/ec/):
  reaction_enzyme_cost.tsv   per-rxn EC, kcat, MW(kDa), cost, per-cond share
  ec_growth_curve.tsv        biomass frac-of-baseline vs budget B, per condition
  ec_flux_matrix.tsv         rxn x condition proteome-based flux at operating B
  ec_saturated_enzymes.tsv   budget-limiting reactions per condition
  ec_report.md               narrative + numbers
"""
import sys, json, re
from pathlib import Path
import numpy as np
import pandas as pd
import cobra

ROOT = Path("/home/janakae/fsp237")
sys.path.insert(0, str(ROOT / "rna_seq_integration" / "src"))
sys.path.insert(0, str(ROOT / "simulations"))
from gpr_expression import aggregate, parse_gpr          # noqa: E402
from run_simulation_panel import apply_media              # noqa: E402
from phase1_eflux import PDB, CONDITIONS                  # reuse media + condition map

MODEL = ROOT / "simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json"
MEANS = ROOT / "proteomics_integration/outputs/proteome_condition_means.tsv"
GENE_EC = ROOT / "proteomics_integration/outputs/gene_ko_ec.tsv"
FASTA = ROOT / "rna_seq_integration/infection_stage_flux_analysis/data/fsp237_proteome.faa"
OUT = ROOT / "proteomics_integration/outputs/ec"
OUT.mkdir(parents=True, exist_ok=True)
BIO = "bio_gsm"
FLUX_EPS = 1e-6

# average residue masses (Da)
AA = dict(A=71.0788,R=156.1875,N=114.1038,D=115.0886,C=103.1388,E=129.1155,
          Q=128.1307,G=57.0519,H=137.1411,I=113.1594,L=113.1594,K=128.1741,
          M=131.1926,F=147.1766,P=97.1167,S=87.0782,T=101.1051,W=186.2132,
          Y=163.1760,V=99.1326)
WATER = 18.01524

# EC-class kcat prior (s^-1); Bar-Even 2011 median ~10, refined by first digit
KCAT_BY_CLASS = {"1":13.7, "2":28.0, "3":42.0, "4":15.0, "5":12.0, "6":15.0, "7":20.0}
KCAT_DEFAULT = 10.0


def gene_mw():
    """Da per gene from sequence (sum residue masses + water)."""
    mw = {}
    gid, seq = None, []
    def flush():
        if gid and seq:
            s = "".join(seq)
            mw[gid] = sum(AA.get(c, 110.0) for c in s) + WATER
    for line in FASTA.open():
        if line.startswith(">"):
            flush(); gid = line[1:].split()[0]; seq = []
        else:
            seq.append(line.strip())
    flush()
    return mw


def main():
    model = cobra.io.load_json_model(str(MODEL))
    means = pd.read_csv(MEANS, sep="\t").set_index("gene")
    ec_tab = pd.read_csv(GENE_EC, sep="\t").set_index("gene")
    gene_ec = {g: str(e).split(";")[0] for g, e in ec_tab["ec"].dropna().items() if str(e)}
    mw = gene_mw()
    med_mw = float(np.median(list(mw.values())))
    print(f"MW: {len(mw)} genes (median {med_mw/1000:.1f} kDa); EC: {len(gene_ec)} genes")

    def rxn_kcat_ec(gpr):
        ecs = [gene_ec[g] for g in {t for cl in parse_gpr(gpr) for t in cl} if g in gene_ec]
        if not ecs:
            return KCAT_DEFAULT, ""
        # representative EC = most common; kcat = mean of class priors
        from collections import Counter
        rep = Counter(ecs).most_common(1)[0][0]
        kcats = [KCAT_BY_CLASS.get(e.split(".")[0], KCAT_DEFAULT) for e in ecs]
        return float(np.mean(kcats)), rep

    def rxn_mw(gpr):
        """min over OR-clauses of sum over AND-tokens of gene MW (Da)."""
        clauses = parse_gpr(gpr)
        if not clauses:
            return None
        clause_mw = []
        for cl in clauses:
            s = sum(mw.get(g, med_mw) for g in cl)
            clause_mw.append(s)
        return min(clause_mw) if clause_mw else None

    # ---- per-reaction cost + per-condition abundance share ----
    lin = {lab: {g: float(2.0**v) for g, v in means[col].dropna().items()}
           for lab, col, _ in CONDITIONS}
    rows = []
    for r in model.reactions:
        gpr = r.gene_reaction_rule or ""
        if not gpr.strip():
            continue
        m_da = rxn_mw(gpr)
        if m_da is None or m_da <= 0:
            continue
        kcat, rep_ec = rxn_kcat_ec(gpr)
        cost = m_da / kcat                        # g protein per unit flux (relative)
        rec = dict(rxn_id=r.id, name=r.name or "", subsystem=r.subsystem or "",
                   ec=rep_ec, kcat=round(kcat, 2), mw_kDa=round(m_da/1000, 2),
                   cost=cost, reversible=(r.lower_bound < 0))
        for lab, _, _ in CONDITIONS:
            ab, nt, ne = aggregate(gpr, lin[lab])
            rec[f"abund_{lab}"] = ab if (ne > 0 and np.isfinite(ab)) else np.nan
        rows.append(rec)
    cdf = pd.DataFrame(rows).set_index("rxn_id")

    # measured reactions = those with abundance in a given condition
    def shares(lab):
        a = cdf[f"abund_{lab}"].dropna()
        return (a / a.sum()).to_dict()

    cost = cdf["cost"].to_dict()
    kcat_d = cdf["kcat"].to_dict()
    mw_d = (cdf["mw_kDa"] * 1000).to_dict()

    # capacity per unit budget: cap_i(cond) = share_i / cost_i  = share*kcat/MW
    cap = {}
    for lab, _, _ in CONDITIONS:
        sh = shares(lab)
        cap[lab] = {rid: sh[rid] / cost[rid] for rid in sh if cost.get(rid, 0) > 0}
    cdf.reset_index().to_csv(OUT / "reaction_enzyme_cost.tsv", sep="\t", index=False)
    print(f"catalyzed reactions costed: {len(cdf)}")

    # ---- budget sweep ----
    # pick B range so ceilings span sub- to super-limiting. cap ~ share/cost;
    # median cap*B should reach ~ typical flux (O(1)) near the knee.
    med_cap = np.median([v for d in cap.values() for v in d.values()])
    print(f"median cap (per unit B): {med_cap:.3e}")
    B_grid = [1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8, 3e8, 1e9, 3e9, 1e10]
    baseline = {}
    curve_rows = []
    for lab, col, media_key in CONDITIONS:
        with model:
            apply_media(model, PDB[media_key], aerobic=True)
            model.objective = BIO
            baseline[lab] = float(model.optimize().objective_value or 0.0)
        for B in B_grid:
            with model:
                apply_media(model, PDB[media_key], aerobic=True)
                model.objective = BIO
                for rid, capv in cap[lab].items():
                    r = model.reactions.get_by_id(rid)
                    ceil = B * capv
                    if r.upper_bound > 0:
                        r.upper_bound = min(r.upper_bound, ceil)
                    if r.lower_bound < 0:
                        r.lower_bound = max(r.lower_bound, -ceil)
                sol = model.optimize()
                bio = float(sol.objective_value or 0.0) if sol.status == "optimal" else 0.0
            curve_rows.append(dict(condition=lab, budget=B, biomass=round(bio, 6),
                frac_of_baseline=round(bio/baseline[lab], 4) if baseline[lab] else None))
    curve = pd.DataFrame(curve_rows)
    curve.to_csv(OUT / "ec_growth_curve.tsv", sep="\t", index=False)

    # each condition evaluated at ITS OWN knee: smallest B reaching >=50% baseline.
    # (a single global B leaves the low-demand conditions unconstrained, hiding
    # their binding enzymes; the per-condition knee is where each is enzyme-limited.)
    piv = curve.pivot(index="budget", columns="condition", values="frac_of_baseline")
    B_knee = {}
    for lab, _, _ in CONDITIONS:
        s = piv[lab]
        hit = s[s >= 0.5]
        B_knee[lab] = float(hit.index.min()) if len(hit) else float(s.index.max())
    print("per-condition knee budgets:", {k: f"{v:.1e}" for k, v in B_knee.items()})

    # ---- per-condition flux + capacity utilization at each condition's knee ----
    flux_mat = pd.DataFrame(index=[r.id for r in model.reactions])
    flux_mat["name"] = [r.name or "" for r in model.reactions]
    flux_mat["subsystem"] = [r.subsystem or "" for r in model.reactions]
    util_rows = []
    for lab, col, media_key in CONDITIONS:
        B = B_knee[lab]
        with model:
            apply_media(model, PDB[media_key], aerobic=True)
            model.objective = BIO
            ceilings = {}
            for rid, capv in cap[lab].items():
                r = model.reactions.get_by_id(rid)
                ceil = B * capv
                ceilings[rid] = ceil
                if r.upper_bound > 0:
                    r.upper_bound = min(r.upper_bound, ceil)
                if r.lower_bound < 0:
                    r.lower_bound = max(r.lower_bound, -ceil)
            sol = cobra.flux_analysis.pfba(model) if model.optimize().status == "optimal" else model.optimize()
        flux_mat[f"flux_{lab}"] = [float(sol.fluxes.get(r.id, 0.0)) for r in model.reactions]
        # capacity utilization = |flux| / ceiling for every catalyzed reaction
        for rid, ceil in ceilings.items():
            f = float(sol.fluxes.get(rid, 0.0))
            if ceil <= FLUX_EPS:
                continue
            util = abs(f) / ceil
            if util < 0.01:
                continue
            util_rows.append(dict(condition=lab, rxn_id=rid,
                name=cdf.at[rid, "name"], subsystem=cdf.at[rid, "subsystem"],
                ec=cdf.at[rid, "ec"], kcat=cdf.at[rid, "kcat"],
                mw_kDa=cdf.at[rid, "mw_kDa"], flux=round(f, 5),
                ceiling=round(ceil, 5), utilization=round(util, 4),
                saturated=bool(util >= 0.98)))
    flux_mat.to_csv(OUT / "ec_flux_matrix.tsv", sep="\t")
    sat = pd.DataFrame(util_rows).sort_values(
        ["condition", "utilization"], ascending=[True, False]) if util_rows else pd.DataFrame()
    sat.to_csv(OUT / "ec_saturated_enzymes.tsv", sep="\t", index=False)

    # ---- report ----
    n_sat = {lab: int(((sat.condition == lab) & sat.saturated).sum()) if not sat.empty else 0
             for lab, _, _ in CONDITIONS}
    rep = [
        "# Phase 2 — proteome-allocated enzyme-capacity model\n",
        f"Model {MODEL.name}. Catalyzed reactions costed: **{len(cdf)}**. "
        f"kcat = EC-class prior (Bar-Even 2011; median {KCAT_DEFAULT} s⁻¹, refined "
        "by EC first digit) — order-of-magnitude, swappable for BRENDA/DLKcat. "
        f"MW from FSP237 proteome sequences (median {med_mw/1000:.1f} kDa).\n",
        "## Enzyme-budget growth curve (frac of own-medium baseline)\n",
        "```\n" + piv.sort_index(ascending=False).to_string() + "\n```\n",
        "Each condition reaches baseline at a different budget — richer medium "
        "supports more growth but demands more enzyme. Per-condition knee (≥50% "
        f"baseline): PDA {B_knee['PDA']:.1e}, half {B_knee['half']:.1e}, "
        f"onetenth {B_knee['onetenth']:.1e}.\n",
        "## Budget-limiting enzymes at each condition's knee\n",
        f"Saturated (utilization ≥0.98) — PDA: {n_sat['PDA']}, half: {n_sat['half']}, "
        f"onetenth: {n_sat['onetenth']}\n",
    ]
    if not sat.empty:
        top = sat.sort_values("utilization", ascending=False).head(15)
        rep.append("Top enzymes by capacity utilization:\n```\n" +
                   top[["condition","rxn_id","name","ec","kcat","flux","utilization"]].to_string(index=False)
                   + "\n```\n")
    rep.append("## Outputs\n- `reaction_enzyme_cost.tsv`, `ec_growth_curve.tsv`, "
               "`ec_flux_matrix.tsv`, `ec_saturated_enzymes.tsv`")
    (OUT / "ec_report.md").write_text("\n".join(rep))
    print("saturated counts:", n_sat)
    print("wrote outputs to", OUT)


if __name__ == "__main__":
    main()
