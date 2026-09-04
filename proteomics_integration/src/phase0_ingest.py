#!/usr/bin/env python3
"""Phase 0 — ingest, QC, harmonize the FSP237 MPLEx proteomics for FBA integration.

Inputs (read-only):
  proteomics-fsp237/Cs_MPLEx_proteins_mycoCosum/
    protein_wide_filtered_imputed.csv   (values used for modeling; 0 NaN)
    protein_wide_filtered.csv           (real per-condition detection; has NaN)
    {half_PDA_vs_PDA,OneTenth_PDA_vs_PDA,OneTenth_PDA_vs_half_PDA}.csv  (limma DE)
    annota_final3_new.xlsx              (gene_ -> uniprot/ko/pfam/desc)
  simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version10_*.json  (V10 GPRs)

Outputs (proteomics_integration/outputs/):
  proteome_condition_means.tsv   gene x {PDA,half,onetenth} log2 means + detection + DE
  proteome_per_replicate.tsv     tidy per-replicate values with clean condition labels
  gene_annotation_map.tsv        gene_ -> annotation + in_model flag
  coverage_report.md             QC + model/annotation coverage

No model, media, or FBA changes. Pure data prep.
"""
import json, textwrap
from pathlib import Path
import numpy as np
import pandas as pd
import cobra

ROOT = Path("/home/janakae/fsp237")
PROT = ROOT / "proteomics-fsp237" / "Cs_MPLEx_proteins_mycoCosum"
OUT  = ROOT / "proteomics_integration" / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
MODEL = ROOT / "simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json"

# canonical condition label -> (sample-name prefix, matched PDB media key)
CONDITIONS = {
    "PDA":      ("fullstrengthPDA", "19_pdb_baseline"),
    "half":     ("halfPDA",         "20_pdb_half"),
    "onetenth": ("onetenthPDA",     "21_pdb_onetenth"),
}

def sample_cols(df, prefix):
    return [c for c in df.columns if c.startswith(prefix)]

# ---------------------------------------------------------------- load matrices
imp = pd.read_csv(PROT / "protein_wide_filtered_imputed.csv").rename(columns={"Unnamed: 0": "gene"})
filt = pd.read_csv(PROT / "protein_wide_filtered.csv").rename(columns={"Unnamed: 0": "gene"})
imp = imp.set_index("gene"); filt = filt.set_index("gene")
assert (imp.index == filt.index).all(), "imputed/filtered gene order differs"
genes = imp.index.tolist()

# ---------------------------------------------------------------- per-rep tidy
long = imp.reset_index().melt(id_vars="gene", var_name="sample", value_name="log2_imputed")
long["condition"] = long["sample"].str.replace(r"_\d+$", "", regex=True).map(
    {p: c for c, (p, _) in CONDITIONS.items()})
long["replicate"] = long["sample"].str.extract(r"_(\d+)$").astype(int)
long["media_key"] = long["condition"].map({c: m for c, (_, m) in CONDITIONS.items()})
# attach real-detection flag from the non-imputed matrix
fl_long = filt.reset_index().melt(id_vars="gene", var_name="sample", value_name="log2_raw")
long = long.merge(fl_long, on=["gene", "sample"], how="left")
long["detected"] = long["log2_raw"].notna()
long = long[["gene", "condition", "media_key", "replicate", "sample",
             "log2_imputed", "log2_raw", "detected"]]
long.to_csv(OUT / "proteome_per_replicate.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- QC (per sample)
qc = []
for cond, (pfx, media) in CONDITIONS.items():
    for c in sample_cols(imp, pfx):
        v = imp[c].values
        rawv = filt[c].values
        qc.append(dict(sample=c, condition=cond,
                       median=np.median(v), mean=np.mean(v),
                       min=np.min(v), max=np.max(v),
                       n_detected=int(np.isfinite(rawv).sum()),
                       n_imputed=int(np.isnan(rawv).sum())))
qc = pd.DataFrame(qc)

# ---------------------------------------------------------------- condition means
means = pd.DataFrame(index=genes)
for cond, (pfx, media) in CONDITIONS.items():
    cols = sample_cols(imp, pfx)
    means[f"{cond}_mean"] = imp[cols].mean(axis=1)
    means[f"{cond}_sd"]   = imp[cols].std(axis=1)
    # genuine detection count per condition (from non-imputed)
    means[f"{cond}_ndet"] = filt[sample_cols(filt, pfx)].notna().sum(axis=1)
means.index.name = "gene"

# ---------------------------------------------------------------- merge DE contrasts
de_files = {
    "half_vs_PDA":      "half_PDA_vs_PDA.csv",
    "onetenth_vs_PDA":  "OneTenth_PDA_vs_PDA.csv",
    "onetenth_vs_half": "OneTenth_PDA_vs_half_PDA.csv",
}
for tag, fn in de_files.items():
    d = pd.read_csv(PROT / fn).rename(columns={"Unnamed: 0": "gene"}).set_index("gene")
    means[f"logFC_{tag}"]   = d["logFC"]
    means[f"adjP_{tag}"]    = d["adj.P.Val"]
means = means.reset_index()

# ---------------------------------------------------------------- annotation map
an = pd.read_excel(PROT / "annota_final3_new.xlsx").rename(columns={
    "ID": "gene",
    "uniprot_chig_id": "uniprot_chig_id",
    "Rank": "annot_rank",
    "ko_id": "ko_id",
    "kegg_hit": "kegg_hit",
    "pfam_hits": "pfam_hits",
    "UP000092177_C_higginsianum description": "chig_description",
    "Identified in proteomics data?": "identified_in_proteomics",
    "unique to one-tenth PDA": "unique_onetenth",
})
an["unique_onetenth"] = an["unique_onetenth"].notna()
# EC is deferred to Phase 3 (kofam on poplar); placeholder column so schema is stable
an["ec_number"] = pd.NA
an["ec_source"] = pd.NA

# ---------------------------------------------------------------- model coverage
m = cobra.io.load_json_model(str(MODEL))
model_genes = {g.id for g in m.genes}
gpr_rxns = [r for r in m.reactions if r.gene_reaction_rule.strip()]

prot_genes = set(genes)
an["in_model"] = an["gene"].isin(model_genes)
an["detected_in_mplex"] = an["gene"].isin(prot_genes)

# keep annotation for all annotated genes, but tag which are in model / detected
annot_map = an[["gene", "uniprot_chig_id", "annot_rank", "ko_id", "kegg_hit",
                "pfam_hits", "ec_number", "ec_source", "chig_description",
                "identified_in_proteomics", "unique_onetenth",
                "in_model", "detected_in_mplex"]]
annot_map.to_csv(OUT / "gene_annotation_map.tsv", sep="\t", index=False)

# attach a couple of annotation cols onto the means table for convenience
means = means.merge(annot_map[["gene", "ko_id", "in_model", "uniprot_chig_id"]],
                    on="gene", how="left")
means.to_csv(OUT / "proteome_condition_means.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- coverage stats
inter_genes = model_genes & prot_genes
gpr_covered = sum(1 for r in gpr_rxns if {g.id for g in r.genes} & prot_genes)
# per-condition: model genes with >=1 real detection in that condition
per_cond_model = {}
for cond, (pfx, media) in CONDITIONS.items():
    det_here = set(filt.index[filt[sample_cols(filt, pfx)].notna().any(axis=1)])
    per_cond_model[cond] = len(model_genes & det_here)

ko_in_model = an[(an.in_model) & an.ko_id.notna()]["gene"].nunique()
detected_model_genes_no_anno = len(inter_genes - set(an.gene))

# ---------------------------------------------------------------- write report
rep = f"""# Phase 0 — proteomics ingest / QC / coverage report

Generated by `proteomics_integration/src/phase0_ingest.py`. Values are TMT **log2**
abundance (relative). Modeling values come from the **imputed** matrix (0 NaN);
genuine per-condition **detection** comes from the non-imputed matrix.

## Condition mapping (proteomics ↔ simulated PDB media)

| Proteomics condition | sample prefix | reps | matched FBA media key |
|---|---|---|---|
| PDA (full strength)  | `fullstrengthPDA` | 5 | `19_pdb_baseline` |
| half-PDA             | `halfPDA`         | 5 | `20_pdb_half` |
| one-tenth-PDA        | `onetenthPDA`     | 5 | `21_pdb_onetenth` |

## Matrix

- Proteins quantified: **{len(genes)}** (native `gene_NNNN` IDs)
- Samples: **{imp.shape[1]}** (3 conditions × 5 reps)
- Imputed matrix NaN cells: **{int(imp.isna().sum().sum())}**  (analysis-ready)
- Non-imputed NaN cells: **{int(filt.isna().sum().sum())}** / {filt.size} = {100*filt.isna().sum().sum()/filt.size:.1f}% (imputed for modeling)

## Per-sample QC (log2)

{qc.round(3).to_markdown(index=False)}

Column medians span {qc['median'].min():.2f}–{qc['median'].max():.2f}
(spread {qc['median'].max()-qc['median'].min():.2f} log2) — {'consistent (no renormalization needed)' if qc['median'].max()-qc['median'].min() < 1.0 else 'CHECK: medians differ >1 log2, consider renormalization'}.

## Model coverage (V10: {len(model_genes)} genes, {len(gpr_rxns)} GPR reactions)

- Model genes detected in MPLEx: **{len(inter_genes)} / {len(model_genes)} = {100*len(inter_genes)/len(model_genes):.1f}%**
- GPR reactions with ≥1 detected protein: **{gpr_covered} / {len(gpr_rxns)} = {100*gpr_covered/len(gpr_rxns):.1f}%**
- Model genes genuinely detected **per condition** (non-imputed):
  - PDA: {per_cond_model['PDA']}   half: {per_cond_model['half']}   onetenth: {per_cond_model['onetenth']}
- Model genes with no row in annotation xlsx: {detected_model_genes_no_anno}

## Annotation coverage

- Annotated genes total: **{len(an)}**  (proteome-wide annotation table)
- With KO (`ko_id`): **{an.ko_id.notna().sum()}**  |  with Pfam: **{an.pfam_hits.notna().sum()}**
- **Model** genes with KO: **{ko_in_model} / {len(model_genes)} = {100*ko_in_model/len(model_genes):.1f}%**
- EC numbers: **deferred to Phase 3** — run kofam (poplar) on the proteome FASTA to
  fill KO for the {len(model_genes)-ko_in_model} model genes lacking KO, then KO→EC.

## Differential expression (limma, adjP<0.05)

| contrast | # significant |
|---|---|
| half vs PDA | {(pd.read_csv(PROT/'half_PDA_vs_PDA.csv')['adj.P.Val']<0.05).sum()} |
| onetenth vs PDA | {(pd.read_csv(PROT/'OneTenth_PDA_vs_PDA.csv')['adj.P.Val']<0.05).sum()} |
| onetenth vs half | {(pd.read_csv(PROT/'OneTenth_PDA_vs_half_PDA.csv')['adj.P.Val']<0.05).sum()} |

## Outputs written

- `outputs/proteome_condition_means.tsv` — {means.shape[0]} genes × {means.shape[1]} cols (means, sd, per-cond detection, DE logFC/adjP, ko_id, in_model)
- `outputs/proteome_per_replicate.tsv` — tidy {long.shape[0]} rows (gene×sample)
- `outputs/gene_annotation_map.tsv` — {annot_map.shape[0]} genes × {annot_map.shape[1]} cols
- `outputs/coverage_report.md` — this file

## Next (Phase 1)

Feed `proteome_condition_means.tsv` (per-condition means) through
`rna_seq_integration/src/gpr_expression.py::aggregate` and
`stage2_gimme_imat_eflux.py` to run E-Flux/GIMME/iMAT per condition on the matched
PDB media, vs pFBA. Expression now varies by condition — the key novelty.
"""
(OUT / "coverage_report.md").write_text(rep)
print(rep)
print("\n[phase0] wrote 4 outputs to", OUT)
