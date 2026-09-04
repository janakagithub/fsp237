#!/opt/env/modelseed/bin/python
"""Assemble atp-safe/proteomics_payload.json from Stage 1 (E-Flux) and Stage 2
(enzyme-capacity) outputs. Pure pandas — no model load. Fetched lazily by the
Proteomics tab, mirroring the infection-stage payload pattern."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/janakae/fsp237")
EF = ROOT / "proteomics_integration/outputs/eflux"
EC = ROOT / "proteomics_integration/outputs/ec"
OUT = ROOT / "atp-safe/proteomics_payload.json"
CONDS = ["PDA", "half", "onetenth"]


def clean(v):
    if isinstance(v, float) and (np.isnan(v) or np.isinf(v)):
        return None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    return v


def recs(df):
    return [{k: clean(v) for k, v in row.items()} for row in df.to_dict("records")]


# ---- Stage 1: E-Flux ----
ef_sum = pd.read_csv(EF / "eflux_summary.tsv", sep="\t")
ef_tit = pd.read_csv(EF / "eflux_titration.tsv", sep="\t")
ef_bot = pd.read_csv(EF / "eflux_bottlenecks.tsv", sep="\t")
ef_conc = json.loads((EF / "eflux_concordance_stats.json").read_text())
ef_fm = pd.read_csv(EF / "eflux_flux_matrix.tsv", sep="\t")

# first bottleneck per condition = the reaction(s) binding at the largest
# capacity where any reaction is throttled
first_bot = {}
for c in CONDS:
    sub = ef_bot[ef_bot.condition == c]
    if len(sub):
        capmax = sub.capacity.max()
        first_bot[c] = sorted(set(sub[sub.capacity == capmax].rxn_id))

# top E-Flux bottleneck rows (dedup by rxn, keep highest-capacity appearance)
bot_top = (ef_bot.sort_values("capacity", ascending=False)
           .drop_duplicates(["condition", "rxn_id"])
           .sort_values(["condition", "capacity"], ascending=[True, False])
           .head(60))[["condition", "capacity", "rxn_id", "name", "subsystem",
                       "flux", "ceiling", "expr_norm"]]

# ---- Stage 2: enzyme capacity ----
ec_curve = pd.read_csv(EC / "ec_growth_curve.tsv", sep="\t")
ec_util = pd.read_csv(EC / "ec_saturated_enzymes.tsv", sep="\t")
ec_cost = pd.read_csv(EC / "reaction_enzyme_cost.tsv", sep="\t")
ec_fm = pd.read_csv(EC / "ec_flux_matrix.tsv", sep="\t")

knees = {}
piv = ec_curve.pivot(index="budget", columns="condition", values="frac_of_baseline")
for c in CONDS:
    hit = piv[c][piv[c] >= 0.5]
    knees[c] = float(hit.index.min()) if len(hit) else float(piv.index.max())

# coverage stats from the enzyme-cost table (one row per catalyzed reaction)
n_costed = len(ec_cost)
n_ec = int(ec_cost["ec"].astype(str).str.strip().replace("nan", "").astype(bool).sum())
n_prot = {c: int(ec_cost[f"abund_{c}"].notna().sum()) for c in CONDS}

# top costly enzymes (MW/kcat) — the enzymes that "cost the most" per unit flux
cost_top = (ec_cost.sort_values("cost", ascending=False).head(40)
            [["rxn_id", "name", "subsystem", "ec", "kcat", "mw_kDa", "cost"]].copy())
cost_top["cost"] = cost_top["cost"].round(1)

# proteome-based flux differentiators across conditions (Stage 2, at knees)
fcols = [f"flux_{c}" for c in CONDS]
fm = ec_fm.copy()
fm[fcols] = fm[fcols].fillna(0.0)
fm["absmax"] = fm[fcols].abs().max(axis=1)
fm["range"] = fm[fcols].max(axis=1) - fm[fcols].min(axis=1)
fm_diff = fm[fm.absmax > 1e-6].sort_values("range", ascending=False).head(45)
fm_diff = fm_diff.rename(columns={ec_fm.columns[0]: "rxn_id"})
flux_diff = fm_diff[["rxn_id", "name", "subsystem"] + fcols].copy()
for c in fcols:
    flux_diff[c] = flux_diff[c].round(4)

payload = {
    "meta": {
        "model": "fsp237_gapfilled_Version10",
        "conditions": CONDS,
        "cond_labels": {"PDA": "full PDB", "half": "½ PDB", "onetenth": "1/10 PDB"},
        "media_keys": {"PDA": "19_pdb_baseline", "half": "20_pdb_half",
                       "onetenth": "21_pdb_onetenth"},
        "n_catalyzed": n_costed,
        "n_ec": n_ec,
        "ec_pct": round(100 * n_ec / n_costed, 1),
        "n_prot_covered": n_prot,
        "n_proteins_dataset": 4981,
        "n_model_genes_covered": 832,
        "n_model_genes": 1274,
        "kcat_source": "EC-class prior (Bar-Even 2011; median 10 s⁻¹, refined by EC first digit)",
    },
    "eflux": {
        "summary": recs(ef_sum),
        "titration": recs(ef_tit),
        "concordance": ef_conc,
        "first_bottlenecks": first_bot,
        "bottlenecks": recs(bot_top),
    },
    "ec": {
        "growth_curve": recs(ec_curve),
        "knees": {k: clean(v) for k, v in knees.items()},
        "utilization": recs(ec_util),
        "cost_top": recs(cost_top),
        "flux_diff": recs(flux_diff),
    },
}

OUT.write_text(json.dumps(payload, separators=(",", ":")))
print(f"wrote {OUT} ({OUT.stat().st_size/1024:.0f} KB)")
print(f"coverage: {n_costed} catalyzed rxns, {n_ec} EC ({payload['meta']['ec_pct']}%)")
print(f"prot-covered: {n_prot}; knees: {knees}")
print(f"eflux bottleneck rows: {len(bot_top)}; ec util rows: {len(ec_util)}; "
      f"cost_top: {len(cost_top)}; flux_diff: {len(flux_diff)}")
