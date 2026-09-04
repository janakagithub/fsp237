#!/opt/env/modelseed/bin/python3
"""Phase 5 — proteome differential-expression & pathway analysis (relative-safe).

Everything here rides the RELIABLE axis of TMT: within-protein fold changes
across the three media (limma logFC / adjP already in proteome_condition_means).
No cross-protein absolute-abundance assumption is used, so these are the
best-supported proteomics readouts. Mirrors the RNA-seq stage-5 analyses and
reuses the shared 23-bucket pathway classifier (pathway_assignment.tsv).

Outputs (proteomics_integration/outputs/pathway/):
  gene_de.tsv                model-gene DE table (3 contrasts, sig flags)
  pathway_fc_matrix.tsv      pathway x contrast mean reaction logFC   [Tier 1.2]
  pathway_enrichment.tsv     Fisher: pathway enriched for DE proteins  [enrichment]
  reporter_metabolites.tsv   Patil-Nielsen on proteome logFC           [Tier 1.4]
  pathway_flux_concordance.tsv  proteome-FC vs GIMME flux-change sign   [Tier 1.3]
  orphan_gpr_gaps.tsv        GPR reactions with no protein detected     [Tier 3.7]
"""
import math
import sys
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, norm

ROOT = Path("/home/janakae/fsp237")
sys.path.insert(0, str(ROOT / "rna_seq_integration/src"))
from gpr_expression import parse_gpr                        # noqa: E402

MODEL_PATH = (ROOT / "simulations/gapfill_v1_v2/models/"
              "fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json")
MEANS = ROOT / "proteomics_integration/outputs/proteome_condition_means.tsv"
PW_ASSIGN = ROOT / "rna_seq_integration/outputs/pathway_assignment.tsv"
GIMME_FM = ROOT / "proteomics_integration/outputs/context/gimme_flux_matrix.tsv"
EC_COST = ROOT / "proteomics_integration/outputs/ec/reaction_enzyme_cost.tsv"
OUT = ROOT / "proteomics_integration/outputs/pathway"
OUT.mkdir(parents=True, exist_ok=True)

CONTRASTS = ["half_vs_PDA", "onetenth_vs_PDA", "onetenth_vs_half"]
SIG_Q, SIG_FC = 0.05, 1.0     # DE call: adjP<0.05 & |logFC|>1


def rxn_genes(model):
    return {r.id: {g for cl in parse_gpr(r.gene_reaction_rule or "") for g in cl}
            for r in model.reactions}


def rxn_mean_fc(genes, fc_by_gene):
    """Signed mean logFC over a reaction's GPR genes (relative-safe)."""
    vals = [fc_by_gene[g] for g in genes if g in fc_by_gene and not np.isnan(fc_by_gene[g])]
    return float(np.mean(vals)) if vals else np.nan


def main():
    model = cobra.io.load_json_model(str(MODEL_PATH))
    means = pd.read_csv(MEANS, sep="\t").set_index("gene")
    pw = pd.read_csv(PW_ASSIGN, sep="\t").set_index("rxn_id")["pathway"].to_dict()
    rg = rxn_genes(model)
    model_genes = {g for gs in rg.values() for g in gs}
    print(f"model {len(model.reactions)} rxns / {len(model_genes)} GPR genes; "
          f"proteome {len(means)} genes")

    # ---- 1. gene-level DE (model genes) ----
    de = means[means.index.isin(model_genes)].copy()
    keep = ["ko_id", "uniprot_chig_id"]
    for c in CONTRASTS:
        keep += [f"logFC_{c}", f"adjP_{c}"]
    de = de[[k for k in keep if k in de.columns]].copy()
    for c in CONTRASTS:
        lf, ap = f"logFC_{c}", f"adjP_{c}"
        if lf in de and ap in de:
            de[f"sig_{c}"] = ((de[ap] < SIG_Q) & (de[lf].abs() > SIG_FC)).fillna(False)
    de.round(4).to_csv(OUT / "gene_de.tsv", sep="\t")
    nsig = {c: int(de.get(f"sig_{c}", pd.Series(dtype=bool)).sum()) for c in CONTRASTS}
    print("DE sig proteins (model genes):", nsig)

    # ---- 2. pathway x contrast mean reaction logFC ----
    fc_dicts = {c: means[f"logFC_{c}"].dropna().to_dict() for c in CONTRASTS
                if f"logFC_{c}" in means}
    rxn_fc = {c: {rid: rxn_mean_fc(rg[rid], fc_dicts[c]) for rid in rg}
              for c in fc_dicts}
    pw_rows = {}
    for rid, p in ((r.id, pw.get(r.id, "Other / unassigned")) for r in model.reactions):
        pw_rows.setdefault(p, []).append(rid)
    fc_matrix = []
    for p, rids in pw_rows.items():
        row = {"pathway": p, "n_rxn": len(rids)}
        for c in fc_dicts:
            vals = [rxn_fc[c][rid] for rid in rids if not np.isnan(rxn_fc[c][rid])]
            row[f"meanFC_{c}"] = round(float(np.mean(vals)), 4) if vals else None
            row[f"n_meas_{c}"] = len(vals)
        pw_rows.setdefault(p, rids)
        fc_matrix.append(row)
    fc_df = pd.DataFrame(fc_matrix).sort_values(
        "meanFC_onetenth_vs_PDA", key=lambda s: s.abs(), ascending=False, na_position="last")
    fc_df.to_csv(OUT / "pathway_fc_matrix.tsv", sep="\t", index=False)

    # ---- 3. Fisher enrichment: pathway enriched for DE proteins (onetenth_vs_PDA) ----
    de_genes = set(de[de.get("sig_onetenth_vs_PDA", False) == True].index) if \
        "sig_onetenth_vs_PDA" in de else set()
    enrich = []
    for p, rids in pw_rows.items():
        pg = {g for rid in rids for g in rg[rid]} & model_genes
        if len(pg) < 3:
            continue
        a = len(pg & de_genes); b = len(pg) - a
        c_ = len(de_genes - pg); d_ = len(model_genes - pg - de_genes)
        try:
            odr, pval = fisher_exact([[a, b], [c_, d_]], alternative="greater")
        except Exception:
            odr, pval = float("nan"), 1.0
        enrich.append(dict(pathway=p, n_genes=len(pg), n_de=a,
                           frac_de=round(a / len(pg), 3),
                           odds_ratio=round(float(odr), 3), p_fisher=round(float(pval), 5)))
    en_df = pd.DataFrame(enrich).sort_values("p_fisher")
    # BH q
    if len(en_df):
        ps = en_df["p_fisher"].values; m = len(ps)
        order = np.argsort(ps); bh = np.empty(m)
        for rank, idx in enumerate(order, 1):
            bh[idx] = ps[idx] * m / rank
        en_df["bh_q"] = np.round(np.clip(np.minimum.accumulate(bh[order][::-1])[::-1]
                                         if m else bh, 0, 1), 4) if m else 1.0
    en_df.to_csv(OUT / "pathway_enrichment.tsv", sep="\t", index=False)
    print("top enriched-for-DE pathways:",
          list(en_df.head(5)["pathway"]) if len(en_df) else [])

    # ---- 4. reporter metabolites on proteome logFC (onetenth_vs_PDA) ----
    rscore = rxn_fc.get("onetenth_vs_PDA", {})
    scores = np.array([s for s in rscore.values() if s is not None and not math.isnan(s)])
    rep_rows = []
    if len(scores) >= 20:
        mu, sig = float(scores.mean()), float(scores.std(ddof=1) or 1.0)
        for met in model.metabolites:
            zs = [(rscore[r.id] - mu) / sig for r in met.reactions
                  if r.id in rscore and not math.isnan(rscore[r.id])]
            if len(zs) < 3:
                continue
            z = float(np.mean(zs) * math.sqrt(len(zs)))
            rep_rows.append(dict(metabolite_id=met.id, name=met.name or "",
                                 compartment=met.id.rsplit("_", 1)[-1] if "_" in met.id else "",
                                 k_reactions=len(zs), z_reporter=round(z, 3),
                                 p_one_sided=round(1 - float(norm.cdf(abs(z))), 4)))
    rep_df = (pd.DataFrame(rep_rows).sort_values("z_reporter", key=lambda c: c.abs(),
              ascending=False) if rep_rows else pd.DataFrame())
    rep_df.to_csv(OUT / "reporter_metabolites.tsv", sep="\t", index=False)
    print(f"reporter metabolites: {len(rep_df)}")

    # ---- 5. pathway-resolved proteome-FC vs GIMME flux-change concordance ----
    conc_rows = []
    if GIMME_FM.exists():
        gf = pd.read_csv(GIMME_FM, sep="\t", index_col=0)
        if {"flux_PDA", "flux_onetenth"}.issubset(gf.columns):
            dflux = (gf["flux_onetenth"].fillna(0) - gf["flux_PDA"].fillna(0)).to_dict()
            fc_o = rxn_fc.get("onetenth_vs_PDA", {})
            for p, rids in pw_rows.items():
                agree = tot = 0
                for rid in rids:
                    f, d = fc_o.get(rid), dflux.get(rid)
                    if f is None or d is None or np.isnan(f) or abs(d) < 1e-6 or abs(f) < 0.5:
                        continue
                    tot += 1
                    if np.sign(f) == np.sign(d):
                        agree += 1
                if tot >= 3:
                    conc_rows.append(dict(pathway=p, n_compared=tot, n_agree=agree,
                                          concordance=round(agree / tot, 3)))
    conc_df = (pd.DataFrame(conc_rows).sort_values("concordance", ascending=False)
               if conc_rows else pd.DataFrame())
    conc_df.to_csv(OUT / "pathway_flux_concordance.tsv", sep="\t", index=False)
    print(f"pathway flux-concordance rows: {len(conc_df)}")

    # ---- 6. orphan-GPR gaps: model reactions with GPR but no protein detected ----
    cov = None
    if EC_COST.exists():
        ec = pd.read_csv(EC_COST, sep="\t").set_index("rxn_id")
        acols = [c for c in ec.columns if c.startswith("abund_")]
        cov = ec[acols].notna().any(axis=1).to_dict() if acols else None
    detected = set(means.index)
    gap_rows = []
    for r in model.reactions:
        gs = rg[r.id]
        if not gs:
            continue
        has_prot = bool(gs & detected)
        covered = cov.get(r.id, has_prot) if cov is not None else has_prot
        if not covered:
            gap_rows.append(dict(rxn_id=r.id, name=r.name or "",
                                 pathway=pw.get(r.id, "Other / unassigned"),
                                 n_genes=len(gs),
                                 genes=";".join(sorted(gs))[:120],
                                 subsystem=r.subsystem or ""))
    gap_df = pd.DataFrame(gap_rows)
    gap_df.to_csv(OUT / "orphan_gpr_gaps.tsv", sep="\t", index=False)
    print(f"orphan-GPR gaps (GPR present, no protein detected): {len(gap_df)} / "
          f"{sum(1 for r in model.reactions if rg[r.id])} GPR reactions")


if __name__ == "__main__":
    main()
