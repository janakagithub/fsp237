# infection_stage_flux_analysis

Systematic comparison of **vanilla FBA** and **transcriptomics-constrained fluxes**
(E-Flux, GIMME, iMAT) for FSP237 V10 across the 18-condition media panel, grouped
by anthracnose infection stage. Milestone 1 (analysis backend + one website tab).

All work here is **downstream of the frozen `rna_seq_integration` outputs** — no FBA
or MILP is re-run. The model, the media panel, the existing RNA-seq pipeline, and the
existing RNA-seq website tab are all treated as immutable.

## The one constraint that shapes everything

There is a **single S1 transcriptome** (technical replicates G1/G2/G3) representing
the *pathogenic state*, applied uniformly to **every** medium. Expression does **not**
vary by infection stage — only flux (via the medium) does. Consequences:

- **No differential-expression analysis** (it would be identically zero across stages).
- S1 is used as a **fixed expression prior**. The two questions we *can* answer:
  1. **Which integration method best reproduces S1** (concordance ranking), and
  2. **Which infection-stage medium the model's flux matches S1 most closely**
     (stage alignment).
- Uncertainty is a **gene-level bootstrap** (resample model genes, *not* the 3
  technical replicates). **No p-values** on stage differences — n = 1 biological
  transcriptome. Intervals are descriptive gene-sampling variability only.

## Pipeline (run order)

All scripts use `/opt/env/modelseed/bin/python3` (no CPLEX needed — Stage-2 MILP
outputs are reused, not recomputed). Run from `scripts/`:

| # | script | output |
|---|--------|--------|
| 00 | `00_build_kegg_map.py` | `data/rxn_kegg_map.tsv` |
| 01 | `01_build_unified_matrix.py` | `flux_results/unified_flux_matrix.tsv` (+ all-thresholds, wide pivots) |
| 02 | `02_concordance_metrics.py` | `method_comparison/{concordance_by_method_condition,concordance_composite,method_ranking}.tsv` |
| 03 | `03_stage_alignment.py` | `statistics/{stage_alignment,stage_alignment_robustness,best_aligned_medium,vanilla_fba_stage_view,vanilla_fba_stage_summary}.tsv` |
| 04 | `04_bootstrap_uncertainty.py` | `statistics/{stage_alignment_bootstrap.tsv,bootstrap_summary.json}` |
| 05 | `05_pathway_method_matrix.py` | `pathway_analysis/{pathway_method_matrix,kegg_pathway_matrix,stage_differential_pathways,method_discordant_pathways}.tsv` |
| 06 | `06_build_web_payload.py` | `../../atp-safe/infection_stage_payload.json` (+ `web/` copy) — M1 tab |
| 07 | `07_build_viz_payload.py` | `../../atp-safe/infection_stage_viz.json` (+ `web/` copy) — M2 heatmap / stages / agreement / reaction table |
| 08 | `08_build_escher.py` | `../../atp-safe/{infection_flux_by_condition.json,infection_escher.html}` (+ `web/` JSON copy) — M2 dynamic Escher map |

`_common.py` holds shared paths, constants, and the reusable `concordance()` metric
(a generalization of the metric block in `../../src/stage1_overlay.py`).

## Methods (as they will appear in the manuscript)

### Flux solutions
Four per (condition, O₂), all from the frozen pipeline:
- **pFBA** — vanilla parsimonious FBA; **expression-independent baseline**.
- **E-Flux** (Colijn 2009) — reaction bounds scaled by `min(1, expr/P99)`.
- **GIMME** (Becker & Palsson 2008) — LP; biomass ≥ 90 % max; penalise flux through
  reactions with `expr < lo_thr`. Featured threshold: `default`.
- **iMAT** (Shlomi 2008) — MILP; maximise (#high-active + #low-inactive); biomass
  floor 10 %. Featured threshold: `default`.

### Concordance with S1 (per method × condition)
On GPR-carrying reactions, comparing each method's flux to the fixed S1 prior:
`spearman_expr_vs_flux` (|flux| vs aggregated log2(TPM+1)), `hi_expr_recall`,
`flux_precision_hi_med`, `agreement_score`, and **MCC**/**Jaccard** of
active (|flux| > 1e-6) vs expressed (bin ∈ {hi, med}). The first four replicate
`stage1_summary.tsv` exactly for pFBA (asserted as a regression check, max abs diff
< 1e-3). **Composite** = mean of the four rank-normalised components
(spearman, hi_recall, flux_precision, MCC).

> **iMAT caveat.** iMAT optimises an expression-agreement objective, so scoring it on
> expression agreement is partly circular, and it activates ~1.5× more reactions than
> the other methods (inflating recall, deflating precision). Report its top rank with
> that caveat; the expression-independent pFBA baseline is the neutral comparator.

### Pathway taxonomy — how definitions were built
Two-layer, **hybrid**:
1. **Primary — 23-bucket rule classifier** (`../../src/stage5_pathway_analysis.py::
   classify_pathway`), an order-sensitive cascade: `reaction.subsystem` (populated for
   only 55/1622 reactions) → id-prefix → diagnostic-metabolite pattern → compartment →
   `Other` fallback. Covers **100 %** of the 1622 reactions; first hit wins, no
   double-counting. Assignments audited in `../../outputs/pathway_assignment.tsv`.
2. **Attached — KEGG identifiers**, joined from a local ModelSEEDDatabase clone
   (snapshotted read-only into `data/rxn_kegg_map.tsv`): KEGG reaction IDs from
   `Aliases/Unique_ModelSEED_Reaction_Aliases.txt` (Source = KEGG) and KEGG pathways
   from `reactions.tsv` `pathways` (KEGG segment `rnNNNNN (Name)`). Keyed on the
   ModelSEED base id (compartment suffix stripped). Coverage: **714/1622 (44 %)** KEGG
   reaction, **660/1622 (41 %)** KEGG pathway across **83** KEGG pathways, plus EC for
   794 (49 %). The ~337 custom reactions (EX_/bio/frxn/BiGG-style) have no ModelSEED
   base and thus no KEGG id — they keep the classifier bucket only.

### Uncertainty
Gene-level nonparametric bootstrap (B = 500, seed 1234): resample model genes with
replacement, re-aggregate every GPR (`aggregate` from `../../src/gpr_expression.py`;
AND = min, OR = max), re-bin on the resampled gene distribution, and recompute the
stage ranking using the **expression-independent pFBA flux** (so resampling expression
is self-consistent — E-Flux/GIMME/iMAT flux is a *function of* the point-estimate
expression and is out of bootstrap scope for M1). Report descriptive 95 % percentile
intervals and the fraction of resamples in which each stage ranks first.

## Headline M1 results

- **Method ranking (composite, aerobic):** iMAT > GIMME > pFBA > E-Flux (with the iMAT
  circularity caveat above).
- **Stage alignment:** S1 aligns most closely with **biotrophic** media — first in
  72 % of gene-bootstrap resamples, and the top stage under pFBA, E-Flux, and GIMME
  (iMAT favours pre-infection). This reproduces the earlier RNA-seq-tab finding that
  the best-fitting single media were `07_sucrose_low` and `06_glucose_low`.
- **Most stage-differential pathways (flux):** pentose-sugar catabolism, pectin/Ashwell,
  peroxisomal β-oxidation — the substrate-specific niche-switch pathways.
- **Most method-discordant pathways:** GAM/maintenance, β-oxidation, TCA, glycolysis,
  PPP — central carbon, where the integration methods most disagree on activity.

## Website

A **new top-level tab "Infection-Stage Flux"** in `../../atp-safe/index.html` (purely
additive; the existing RNA-seq tab is untouched). `build_atp_safe_site.py` regenerates
only `reactions.json`, so these edits and the payloads survive rebuilds.

**M1 sub-tabs** (lazy-load `atp-safe/infection_stage_payload.json`): *Overview & Methods*,
*Method Comparison*, *Stage Alignment*, *Pathways*.

**M2 sub-tabs** (lazy-load a second, heavier `atp-safe/infection_stage_viz.json` fetched
only when an M2 sub-tab is first opened, so the M1 tab stays fast):
- *Pathway Heatmap* — pathway × medium active-rate heatmap, method + 23-bucket↔KEGG
  taxonomy toggles (`interpHex` blue ramp). Built by `07_build_viz_payload.py`.
- *Infection Stages* — media-by-stage cards + a stage × method composite-alignment matrix
  + best-aligned-medium table.
- *Flux–Expression Agreement* — stacked agreement-category bars per medium/method over the
  six-category glossary (discordant subset = PRIMED_NOT_USED / ORPHAN_FLUX /
  CONFLICT_FLUX_NO_EXPR).
- *Reaction & Gene Table* — DataTable over all 1622 reactions (pathway, KEGG, EC, GPR,
  log2TPM, per-method active-fraction / mean|flux| / dominant agreement / discordant rate).
- *Escher Map* — ONE dynamic Escher builder page (`infection_escher.html`, cloned from the
  frozen `map_aerobic.html` with a data-driven recolor swapped in) embedded in an iframe;
  a method button-group + medium `<select>` repaint it live from
  `atp-safe/infection_flux_by_condition.json` (built by `08_build_escher.py`, nonzero
  aerobic fluxes only) via `?method=&cond=` URL params. Edge colour + width = |flux|
  (log-scaled); grey = no flux. The frozen `map_aerobic/anaerobic/expression` maps are
  not modified.

## Deferred to later milestones

- **M3** — deep-research biological interpretation (tiered evidence + literature) of the
  top stage-differential / discordant pathways → `deep_research/` + a new sub-tab.
