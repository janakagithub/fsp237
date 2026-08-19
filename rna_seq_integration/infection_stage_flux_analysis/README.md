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
| 09 | `09_build_interp_payload.py` | `../../atp-safe/infection_interp.json` (+ `web/` copy) — M3 biological-interpretation cards (renders `deep_research/pathways/*.md`) |
| 10 | `10_build_gene_annot.py` | `../../atp-safe/infection_gene_annot.json` (+ `web/` copy) — M4 reaction→gene drill-down + pathogenesis layer; also `data/{fsp237_proteome.faa,gene_annotation.tsv,secretome_predictions.tsv}` |

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

**M3 sub-tab** (lazy-loads `atp-safe/infection_interp.json`):
- *Biological Interpretation* — deep-research reading of the six priority pathways (most
  stage-differential + most method-discordant), as filterable cards with a category badge,
  a calibrated confidence badge, a one-liner, an expandable full interpretation, and
  verified key references. Source prose lives in `deep_research/pathways/*.md` (one file per
  pathway), authored with **tiered evidence** (Tier 1 *Colletotrichum*-specific → Tier 2
  related phytopathogens → Tier 3 general fungal/biochemical) per `deep_research/
  RESEARCH_BRIEF.md`. All citations were verified against live NCBI/Europe PMC/Crossref
  records (no fabricated DOIs); `09_build_interp_payload.py` renders the markdown to HTML.

### M3 interpretation highlights (with honest caveats)
- **Plant cell-wall degradation** (pectin/D-galacturonate + pentoses) switches on *only* in
  necrotrophic media — a substrate-driven signature of the tissue-maceration phase, backed by
  *C. sublineola*-specific dual-RNA-seq (Vela 2024) and comparative genomics (Buiate 2017).
  **Annotation flag:** the model's *"Ashwell"* label is the *bacterial oxidative* uronate
  route; fungi use the *reductive* gaaA–gaaD pathway — flux direction holds, but the reaction
  annotation should be checked before publication.
- **Peroxisomal β-oxidation** carries flux *only* pre-infection — matches the storage-lipid →
  β-oxidation → glyoxylate-cycle → appressorial-turgor program (direct *Colletotrichum*
  PEX-mutant + *Magnaporthe* evidence). Direction trustworthy; magnitude not (2nd-most
  method-discordant, single static transcriptome).
- **GAM / maintenance** (top method-discordance, range 1.0) is a **modeling artifact, not
  biology**: the maintenance reaction has no GPR, so iMAT leaves it on (~0.001) in all 18/18
  conditions while pFBA/E-Flux/GIMME force it off — invariant to stage.
- **Central carbon** (glycolysis/TCA/PPP) is uniformly active but the top discordance hotspot
  by network topology + method design (incl. the iMAT active-count inflation caveat).
- **Storage (glycogen/trehalose)** and **fungal cell-wall polysaccharide** give small,
  medium-driven signals that *corroborate* known appressorial-turgor / wall-masking biology
  rather than demonstrate stage-specific regulation (Medium / Medium-Low confidence).

## Milestone 4 — Deep enhancements (genes/pathogenesis, glossary, heatmap-by-stage, brighter Escher)

Four publication-review upgrades, all additive; the RNA-seq tab and the frozen model/panel
stay untouched, and **no FBA/MILP is re-run**.

### Reaction→gene drill-down + pathogenesis layer (`10_build_gene_annot.py`)
- Extracts the proteome via `blastdbcmd -db gpr-update/blast_db/fsp237 -entry all -outfmt '%f'`
  (14,857 proteins; titles carry `gene_NNNN len=… func=<description>`) → `data/fsp237_proteome.faa`
  and `data/gene_annotation.tsv`.
- Loads per-gene S1 expression (`expression-data/S1_normalized_expression.xlsx`, 13,047 genes),
  bins on the **model-gene** quartiles (hi ≥ 6.457 = Q75, lo ≥ 3.121 = Q25 of Mean log2(TPM+1)),
  and maps GPR genes → reactions from `outputs/reaction_expression.tsv`.
- **Two honest layers** (the GEM's GPR genes are metabolic *enzymes*; canonical fungal effectors
  are small secreted proteins that sit *outside* a metabolic gene set — so they are reported
  separately):
  - `virulence_metabolic` — 62 GPR genes in curated virulence enzyme families (cutinase/lipase,
    peroxisomal β-oxidation, ROS detox, chitin/glucan synthesis, P450, melanin, glyoxylate,
    trehalose, pectin/cell-wall degradation).
  - `candidate_effectors` — top 250 **non-metabolic**, highly-expressed genes predicted
    secreted/effector-like (e.g. cerato-platanin, HR-inducing proteins, small hypotheticals) —
    the "highly expressed but not metabolic" layer, explicitly badged as outside the model.
  - `by_reaction` — per-GPR-reaction gene list (gene, func, log2TPM, bin, secreted, effector,
    virulence family) powering the client-side drill-down.
- **Secretome method (documented caveat):** SignalP 6.0 / EffectorP 3.0 are **not runnable here**
  (no EMBOSS pepstats, no torch, no academic license). A transparent **sequence-feature
  heuristic** is used instead — N-terminal Kyte-Doolittle hydrophobic signal-peptide core + small
  + cysteine-rich + non-metabolic function — and flags are labeled *indicative, not definitive*
  in the payload `meta` and on the page. Manuscript follow-up: rerun with SignalP/EffectorP.

### Website (all in `atp-safe/index.html`, additive)
- **Reaction & Gene table:** clicking a reaction id expands a DataTables **child row** listing its
  GPR genes with expression bin + pathogenesis flags (`ensureInfGene` lazy-loads
  `infection_gene_annot.json`).
- **New 11th sub-tab "Genes & Pathogenesis":** metabolic virulence-factor table + non-metabolic
  candidate-effector DataTable + prominent honest caveat callout.
- **Pathway Heatmap** x-axis now grouped by infection stage (pre-infection → biotrophic →
  necrotrophic → cocktail) with `colspan` stage-group headers.
- **Metric glossary + tooltips** (`<abbr>` + collapsible `<details>`) on Point, Bootstrap mean,
  95% CI, Spearman ρ, hi_recall, flux_prec, MCC, composite, agree, active fluxes, Biomass.
- **Escher map** recolored with a brighter multi-stop spectrum (blue→cyan→green→yellow→orange→red)
  with a raised low-flux brightness floor (`08_build_escher.py`), so faint-but-active edges are
  visible; legend gradient updated to match.

## Deferred / follow-ups

- Verify the D-galacturonate reaction annotation ("Ashwell" oxidative vs fungal reductive
  gaaA–gaaD) in the model before manuscript submission.
