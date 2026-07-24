# FSP237 V10 × S1 Transcriptomics Integration — Methods

## Model & data

- **Model**: FSP237 V10 gap-filled + VLCFA-complete + gene-integrated GEM,
  1622 reactions / 1268 metabolites / 1274 genes (frozen for this analysis).
- **Expression**: S1 biological group, 3 technical replicates (G1, G2, G3),
  reported as `log2(TPM + 1)` and reconstructed `TPM`. Gene IDs match the
  model's `gene_NNNN` namespace directly (no ID mapping).
- **Coverage**: 1174 / 1274 model genes (92.2 %) have measured expression;
  the remaining 100 are legacy identifiers (27 `CH63R_*`, 73 yeast `Y*` +
  RefSeq `NP_*` + `SPONT`) that survived the earlier GPR overhaul.

## Stage 0 — gene → reaction expression aggregation

Each reaction's Gene-Protein-Reaction (GPR) rule is parsed into disjunctive
normal form and reduced to a single scalar via
`score(rxn) = max_over_OR( min_over_AND( expr[gene] ) )` — isoenzymes
contribute the loudest expression, complex subunits the quietest. Missing
genes are excluded from the min; a reaction whose entire GPR is missing is
tagged **absent**.

Reactions are placed into four bins on the distribution of `Mean log2(TPM+1)`
restricted to model genes:

| Bin | Definition | Cutoff (this run) |
|---|---|---|
| `hi` | score ≥ 75%ile | ≥ 6.457 |
| `med` | 25%ile ≤ score < 75%ile | [3.121, 6.457) |
| `lo` | 0 < score < 25%ile | (0, 3.121) |
| `absent` | score = 0 or gene missing | — |

## Stage 1 — pFBA × expression overlay

For every (condition, O₂) pair in the 18-condition panel, parsimonious FBA
(pFBA) is solved with the exact medium used by `run_simulation_panel.py`.
Reaction fluxes are joined with the Stage-0 scores and each reaction is
classified:

| Category | Rule |
|---|---|
| `SUPPORTED` | `|v| > ε` AND bin ∈ {hi, med} |
| `WEAK_SUPPORT` | `|v| > ε` AND bin = lo |
| `CONFLICT_FLUX_NO_EXPR` | `|v| > ε` AND bin = absent AND has GPR |
| `ORPHAN_FLUX` | `|v| > ε` AND no GPR (feeds Stage 4) |
| `PRIMED_NOT_USED` | `|v| ≤ ε` AND bin = hi |
| `SILENT_OK` | remainder |

Four condition-fit metrics roll up per (condition, O₂):

- **`agreement_score`** = SUPPORTED / (SUPPORTED + CONFLICT_FLUX_NO_EXPR + PRIMED_NOT_USED).
  Biased toward media that engage more machinery; reported for completeness.
- **`spearman_expr_vs_flux`** — Spearman ρ between the aggregated expression
  and `|v_pFBA|` over GPR'd reactions. Least biased by flux density; used as
  the primary condition-fit ranking on the site.
- **`hi_expr_recall`** — fraction of hi-expression reactions with `|v| > ε`.
- **`flux_precision_hi_med`** — fraction of flux-carrying GPR'd reactions
  whose expression is hi or med.

## Stage 2 — three independent context-specific analyses

No single expression-integration method is universally best across systems
(Blazier & Papin 2012; Machado & Herrgård 2014; recent 2020s benchmarks), so
three techniques run per (condition, O₂) and are compared on the site:

### E-Flux (Colijn 2009)
Continuous. For each reaction with a GPR-derived score, upper (and reversibly
lower) bounds are scaled by `min(1, expr / P99(expr))`. pFBA is solved on the
constrained model. No threshold parameter.

### GIMME (Becker & Palsson 2008)
LP. Biomass is fixed to ≥ 90 % of the unconstrained maximum; the objective
then minimizes Σ (lo_thr − expr[r]) · |v_r| over reactions with `expr < lo_thr`.
Reactions retained despite the penalty are the "must-carry" set. Threshold
sweep: `default (lo = 3.12)`, `strict (lo = 0.88)`, `narrow (lo = 4.99)` on
`log2(TPM+1)`.

### iMAT (Shlomi 2008)
MILP. Reactions with `expr ≥ hi_thr` are class H; `expr < lo_thr` are class L;
otherwise M. Binary indicators are introduced for `|v| ≥ ε` (H) and
`|v| ≤ ε` (L). Objective: maximize `Σ y_H + Σ z_L` under a biomass floor
(≥ 10 % of max) and the standard stoichiometric constraints. Solved with
CPLEX, `mipgap = 0.05`, `timelimit = 60 s`. Two thresholds:
`default` and `strict`. The `narrow` threshold produced a MILP too large for
the 60-s budget and was dropped for iMAT (kept for GIMME).

## Stage 3 — MADE (deferred)

MADE requires ≥ 2 biological groups with differential expression statistics.
The current dataset is single-condition; MADE is bolted on when additional
groups (S2, S3, …) arrive.

## Stage 4 — orphan-reaction expression triage

Reactions flagged `ORPHAN_FLUX` in Stage 1 are grouped across conditions and
ranked by cross-condition flux magnitude. For each, up to three candidate
gene IDs are proposed from within the model: genes appearing in reactions in
the same subsystem, ordered by their own expression score. This does **not**
run BLAST; it produces a prioritised triage list that feeds the future
gene-assignment effort.

## Software

- Python 3.10 / 3.12; cobra 0.29–0.30, pandas 2.2, scipy, openpyxl.
- Solvers: GLPK (LP, default) and CPLEX (MILP, required for tractable iMAT).
- Environments used (this box): `/opt/env/modelseed` for Stages 0, 1, 4;
  `/opt/env/modelseed_cplex` for Stage 2.

## References (methods)

- Becker SA, Palsson BØ. *PLoS Comput Biol* 2008 4:e1000082 — GIMME.
- Shlomi T et al. *Nat Biotechnol* 2008 26:1003 — iMAT.
- Colijn C et al. *PLoS Comput Biol* 2009 5:e1000489 — E-Flux.
- Jerby L, Shlomi T, Ruppin E. *Mol Syst Biol* 2010 — flux–omics integration overview.
- Blazier AS, Papin JA. *Front Physiol* 2012 3:299 — integration-method comparison.
- Machado D, Herrgård M. *PLoS Comput Biol* 2014 10:e1003580 — systematic evaluation.
