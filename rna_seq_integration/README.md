# rna_seq_integration

Transcriptomics × FSP237 V10 GEM integration pipeline. Produces the RNA-seq tab
on the published site and feeds a candidate-gene triage list back into the
gap-fill curation effort.

## Frozen model baseline

- Model: `simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json`
  (1622 rxns / 1268 mets / 1274 genes)
- Panel: the 18-condition medium set from `simulations/run_simulation_panel.py`
- Both are treated as **immutable** for this analysis — do not edit the model
  or the panel from inside this folder.

## Expression dataset

- Path: `/home/janakae/fsp237/expression-data/S1_normalized_expression.xlsx`
  (also reachable via `data/S1_normalized_expression.xlsx` symlink here).
- One biological group **S1** with 3 technical replicates (G1, G2, G3).
- Gene IDs use the same `gene_NNNN` convention as the V10 model, so **no ID
  mapping is required**.
- 92 % of model genes have expression (1174 / 1274). The remaining 100 are
  legacy IDs (27 `CH63R_*`, 73 yeast `Y*` + RefSeq `NP_*` + `SPONT`) that
  survived the earlier GPR overhaul; they surface in the coverage tile so
  they can be cleaned up later.

## Stages

| # | script | tools | outputs |
|---|--------|-------|---------|
| 0 | `src/gpr_expression.py` | GPR parser (reused from `gpr-update/`), pandas | `outputs/reaction_expression.tsv`, `outputs/coverage_summary.json` |
| 1 | `src/stage1_overlay.py` | cobra pFBA + scipy Spearman ρ | `outputs/stage1_overlay/<cond>_<o2>.tsv`, `outputs/stage1_summary.tsv` |
| 2 | `src/stage2_gimme_imat_eflux.py` | cobra + CPLEX (LP + MILP) | `outputs/stage2_eflux/`, `outputs/stage2_gimme/`, `outputs/stage2_imat/`, `outputs/stage2_summary.tsv` |
| 3 | *(deferred)* | MADE / cross-condition DE — awaits S2, S3, … | — |
| 4 | `src/stage4_orphan_validation.py` | plain pandas — reuses Stage-1 outputs | `outputs/stage4_orphan/<cond>_<o2>.tsv`, `outputs/orphan_priority.tsv` |
|   | `src/build_rnaseq_payload.py` | pandas | `outputs/rnaseq_payload.json` (consumed by the site builder) |

### Stage 0 — gene → reaction expression

Aggregation: for each reaction's GPR, `score = max_OR( min_AND( expr[gene] ) )`
using `Mean log2(TPM+1)`. Bins are set on the model-gene distribution:
`hi ≥ 75%ile, med [25%, 75%), lo (0, 25%), absent = 0 or gene missing`.

### Stage 1 — pFBA × expression overlay

Runs pFBA on all 18 × 2 conditions, merges with the Stage-0 scores, and tags
every reaction with an agreement category
(`SUPPORTED / WEAK_SUPPORT / CONFLICT_FLUX_NO_EXPR / PRIMED_NOT_USED /
ORPHAN_FLUX / SILENT_OK`).

Four condition-fit metrics roll up per (condition, O₂):
- `agreement_score` — SUP / (SUP + CFL + PNU)
- `spearman_expr_vs_flux` — Spearman ρ between agg expression and |pFBA flux|
  on GPR'd reactions (**headline metric**; least biased by flux-heavy media)
- `hi_expr_recall` — fraction of hi-expression reactions carrying flux
- `flux_precision_hi_med` — fraction of flux-carrying GPR'd reactions with
  hi/med expression

### Stage 2 — three independent context-specific analyses

Following the manuscript-safe advice that **no single expression-integration
method is universally best**, three techniques run per condition:

- **E-Flux** (Colijn 2009) — continuous. Bounds scaled by
  `min(1, expr / P99(expr))`. Reactions without GPR keep default bounds.
- **GIMME** (Becker & Palsson 2008) — LP. Enforce biomass ≥ 90 % max,
  minimize sum of |flux| weighted by `(lo_thr - expr)` for reactions with
  `expr < lo_thr`. Three thresholds swept (`default / strict / narrow`).
- **iMAT** (Shlomi 2008) — MILP. Reactions binned H / M / L. Maximize
  `#H active + #L inactive`. Two thresholds swept (`default / strict`);
  the `narrow` threshold produces a MILP too large for a 60-s time budget.
  CPLEX with `mipgap=0.05, timelimit=60s` per solve.

### Stage 4 — orphan reactions carrying flux

Filters Stage-1 outputs to `ORPHAN_FLUX` reactions, collects the top-3
subsystem-proximal genes as candidates, ranks by cross-condition flux
magnitude. This does **not** run BLAST — it produces a prioritised list for
the future gene-assignment effort.

## Running the pipeline

```
cd /home/janakae/fungalTemplate/imm904CobraModel/rna_seq_integration

# Stage 0 (safe re-run; no cobra needed)
/opt/env/modelseed/bin/python3 src/gpr_expression.py

# Stage 1 (pFBA — few minutes)
/opt/env/modelseed/bin/python3 src/stage1_overlay.py

# Stage 2 (LP + MILP — 20-40 min with CPLEX)
/opt/env/modelseed_cplex/bin/python3 src/stage2_gimme_imat_eflux.py

# Stage 4 (fast — reuses Stage 1 outputs)
/opt/env/modelseed/bin/python3 src/stage4_orphan_validation.py

# Bundle for the site
/opt/env/modelseed/bin/python3 src/build_rnaseq_payload.py

# Rebuild site
/opt/env/modelseed/bin/python3 /home/janakae/fungalTemplate/imm904CobraModel/fsp237_biomass_extension/build_atp_safe_site.py
```

## Interpreting the site's RNA-seq tab

- **Coverage tiles** — quick sanity numbers.
- **Condition fit** — sort by `Spearman ρ (aer)` to see which panel medium
  best matches the S1 transcriptome. Under S1 the top rows are
  `07_sucrose_low` (invertase → glucose+fructose) and `06_glucose_low`,
  consistent with the biological expectation that S1 is a sugar-utilization
  transcriptome.
- **Stage 2** — sort by `Method` to compare E-Flux / GIMME / iMAT counts
  across the same condition; disagreement between methods is a manuscript
  discussion point, not a bug.
- **Per-reaction** — filter `Bin = hi` + `Agreement = PRIMED_NOT_USED` to see
  highly-expressed reactions whose flux the model does not use — candidates
  for pathway / medium refinement.
- **Orphan reactions** — top rows here feed the next BLAST batch.

## Deferred work

- **MADE / DE** — bolt on when S2 / S3 arrive.
- **10 remaining unassigned gap-fill reaction gene assignments** — will
  consume `outputs/orphan_priority.tsv` as input.
- **Panel extension** (dhurrin, ferulate cross-links) — separate effort;
  the current 18-condition panel is frozen for this analysis.
