# FSP237 × S1 — Pathway-Level Analysis

**Model**: FSP237 V10 GEM (1622 reactions, 1268 metabolites, 1274 genes).
**Expression**: S1 biological group, 3 replicates, `log2(TPM+1)`.
**Reporting date**: July 2026.

## Executive summary

The S1 transcriptome exhibits a canonical **sugar-utilization** metabolic
program: glycolysis, TCA, pentose phosphate, and fatty-acid biosynthesis are
all significantly enriched for highly-expressed reactions (Fisher's exact,
BH-adjusted q ≤ 0.004). This matches the hypothesis that S1 was recovered
from a glucose-minimal-medium-like state — an independent, orthogonal
confirmation of the Spearman-ρ ranking already surfaced on the RNA-seq tab
(top-2 = `07_sucrose_low` and `06_glucose_low`).

Reporter-metabolite analysis (Patil-Nielsen 2005) surfaces mitochondrial
protons, NAD/NADH, thioredoxin, glyceraldehyde-3-phosphate, and malonyl-CoA
as the top nodes whose neighbouring reactions are transcriptionally
activated — again all central-carbon + redox nodes, consistent with active
respiratory metabolism.

## Pathway classification — 20+ heuristic buckets

`reaction.subsystem` is populated for only 55 of 1622 reactions (the
gap-filled specialty pathways). We impose a rule-based classifier
(`src/stage5_pathway_analysis.py::classify_pathway`) that cascades:
subsystem → id-prefix → metabolite pattern → compartment → fallback.
Full per-reaction assignments in `outputs/pathway_assignment.tsv`.

Bucket counts on V10:

| bucket | n reactions | notes |
|---|---:|---|
| Other / unassigned | 576 | mostly generic cofactor + shuttle rxns without a diagnostic metabolite |
| Transport (inter-compartmental) | 256 | `TRP_*`, `tx_*`, or "transporter" in name |
| Exchange | 170 | `EX_*` boundary reactions |
| Amino-acid metabolism | 162 | reactions touching any of the 20 proteinogenic AAs |
| Cofactor / vitamin biosynthesis | 117 | NAD/CoA/FMN/thiamin/folate/biotin/riboflavin/pantothenate/pyridoxal name match |
| Glycolysis / gluconeogenesis | 58 | 15 diagnostic metabolites (G6P, F6P, FBP, DHAP, GAP, PEP, pyruvate…) |
| β-oxidation (peroxisomal) | 54 | subsystem starts "beta-oxidation" OR peroxisomal + acyl-CoA name |
| Fatty-acid biosynthesis | 38 | malonyl-CoA / ACP / long-FA metabolites in non-peroxisomal comps |
| Lipid / membrane | 35 | phosphatidyl-, ergosterol, sphingolipids, cardiolipin, lanosterol |
| TCA cycle | 34 | ≥ 1 TCA-diagnostic metabolite (citrate, α-KG, succinate, fumarate, malate, OAA, succinyl-CoA) |
| Oxidative phosphorylation / ETC | 24 | H⁺ moving between mitochondria + cytosol, or O₂ + H⁺ in mitochondria |
| Nucleotide metabolism | 21 | NTP substrate + nucleobase companion |
| Pentose phosphate pathway | 20 | 6PGL, 6PG, Ru5P, Xu5P, R5P, S7P, E4P |
| Cell-wall polysaccharide | 18 | chitin, α-1,3-glucan, β-1,3-glucan, mannan, UDP-GlcNAc, UDP-glucose |
| Storage (glycogen / trehalose) | 7 | glycogen, trehalose metabolites |
| Melanin / DHN | 7 | subsystem-assigned during biomass extension |
| Biomass | 6 | `bio_gsm`, `bio1`, and the biomass-precursor stoich rxns |
| Pectin / D-galacturonate (Ashwell) | 6 | gap-filled Ashwell pathway |
| Pentose-sugar catabolism | 6 | xylose (XR/XDH/XK) + L-arabinose (Penttilä) pathways |
| Sulfur / nitrogen assimilation | 3 | sulfate/ammonia/nitrate/urea metabolism |
| Compatible solute / mannitol | 2 | MpdA branch |
| GAM / maintenance | 1 | ATP + H₂O → ADP + Pi + H⁺ |
| Sink / demand | 1 | `DM_/SK_` prefix |

The rule set is order-sensitive (first hit wins); no reaction is
double-counted. Alternative-classifier drop-in should be trivial (edit
`classify_pathway`).

## Q1 — Which pathways are active in condition A vs condition B?

The site's **Pathway × condition activity heatmap** (RNA-seq tab) shows,
per pathway × per medium × per O₂ state, the fraction of pathway reactions
carrying flux in pFBA. Enriched-for-hi-expression pathways are stacked at
the top.

Selected shifts between medium classes (mean active-rate across the 3–5
conditions in each class, aerobic):

| pathway | biotrophic (5) | necrotrophic (5) | pre-infection (5) |
|---|---:|---:|---:|
| Glycolysis / gluconeogenesis | high | high | high (as gluconeogenesis from FA / glycerol) |
| TCA cycle | high | high | high |
| β-oxidation (peroxisomal) | ~0 on sugars | ~0 on sugars | **fully on** for palmitate / oleate / hexacosanoate |
| Pentose-sugar catabolism | ~0 | **fully on** for xylose / arabinose | ~0 |
| Pectin / Ashwell | ~0 | **fully on** for D-galacturonate | ~0 |
| Storage (glycogen / trehalose) | ~0 | ~0 | **on** for trehalose |
| Cell-wall polysaccharide | steady across all growable media (biomass-driven) |
| Cofactor / vitamin | steady across all growable media (biosynthesis-driven) |

The heatmap makes the **niche-switch** clear: β-oxidation, Ashwell, and
Penttilä light up only when the corresponding substrate is the sole
carbon source — this validates the direction-locked (degradation-only)
V5/V6 curation and the V1 gap-fill choices, and shows the model does not
routinely leak flux through those pathways.

## Q2 — Pathway enrichment for the S1 transcriptome

Fisher's exact test (one-sided, greater) on `hi vs (med + lo + absent)`
per pathway vs the model background. BH correction across the 23 tested
pathways.

Significant hits (q < 0.05):

| pathway | hi / n | OR | q |
|---|---:|---:|---:|
| **Glycolysis / gluconeogenesis** | 37 / 58 | 5.83 | < 10⁻⁴ |
| **TCA cycle** | 21 / 34 | 5.15 | < 10⁻⁴ |
| **Fatty-acid biosynthesis** | 20 / 38 | 3.52 | 0.001 |
| **Pentose phosphate pathway** | 12 / 20 | 4.69 | 0.004 |
| **Cofactor / vitamin biosynthesis** | 44 / 117 | 1.95 | 0.004 |

Interpretation: **the machinery of aerobic central carbon plus the
biosynthetic charge (fatty-acid, cofactor, plus glycolytic gluconeogenic
demand for cell-wall precursors) is transcriptionally on**. This is the
signature of a well-fed, mycelium-like state — the biological expectation
for a glucose-minimal S1 sample.

Marginal but not enriched (0.05 < q < 0.2): amino-acid metabolism (q =
0.06), cell-wall polysaccharide (q = 0.06). Both are plausibly on but not
statistically distinguishable from the model background.

The full enrichment table lives in `outputs/pathway_enrichment.tsv`; every
tested pathway is also shown on the site's **Pathway detail** table with
its OR, raw p, and BH q.

## Q3 — Reporter metabolites (Patil-Nielsen)

For each metabolite `m` we z-score reaction expression against the
background, average over `m`'s neighbours, and scale by `√k`. Top-10 by
|z| (see `outputs/reporter_metabolites.tsv` for the full list):

| metabolite | comp | k | z | interpretation |
|---|---|---:|---:|---|
| thioredoxin (trdrd / trdox) | c0 | 11 | +6.27 | redox pool; loaded by biosynthesis + PPP + AA-metabolism activity |
| H⁺ (mitochondrial) | m0 | 93 | +5.73 | proton-motive-force node; ETC + TCA activity |
| NADH (mitochondrial) | m0 | 25 | +5.53 | oxidative respiration active |
| NAD⁺ (mitochondrial) | m0 | 26 | +5.37 | same pool, oxidised half |
| glyceraldehyde-3-phosphate | c0 | 8 | +5.03 | glycolytic hub — expression-loaded |
| H⁺ (extracellular) | e0 | 53 | −4.81 | boundary proton — **silent** (S1 is not actively pumping) |
| malonyl-CoA | c0 | 11 | +4.79 | FA-biosynthesis node — matches enrichment result |

The signal is highly consistent with the enrichment result: **redox +
oxidative-respiration + biosynthetic hubs are activated**; extracellular
transport is quiet.

## Q4 — Does the model explain observed growth or phenotype shifts?

Compared pathway-total-flux (Σ|v|) against biomass across the 18 aerobic
conditions using Pearson r. Every biosynthesis / precursor-supply
pathway correlates r ≥ 0.98 with biomass:

| pathway | r | p |
|---|---:|---:|
| Biomass | +1.00 | — |
| Compatible solute / mannitol | +1.00 | < 10⁻⁴ |
| Nucleotide metabolism | +1.00 | < 10⁻⁴ |
| Cell-wall polysaccharide | +0.999 | < 10⁻⁴ |
| Oxidative phosphorylation / ETC | +0.987 | < 10⁻⁴ |

**Manuscript caveat**: these near-perfect correlations are **inherent
FBA collinearity**, not independent evidence. Under pFBA, any biosynthesis
pathway that is on the shortest path to biomass will scale linearly with
biomass. Report this as a sanity check ("the model uses biosynthesis
pathways in proportion to growth demand — no orphan biosynthesis or
unaccounted precursor sinks") rather than as a mechanistic finding.

The **informative** flux-vs-biomass patterns are the pathways whose flux
is **decoupled** from biomass — e.g., β-oxidation (on only for FA
substrates), Ashwell (on only for D-galU), Penttilä (on only for
pentoses). The heatmap surfaces these directly.

## Expression-painted Escher map

New tab **Escher map — expression** (`atp-safe/map_expression.html`) uses
the same iMM904 CCM layout as the flux maps, but paints intensity from
`agg_mean log₂(TPM+1)` instead of |flux|. Palette and legend layout match
the flux maps 1:1 so the three maps can be compared side-by-side (open in
separate tabs).

- Compartment colour (cytosol pink, mito green, ER blue, peroxisome
  orange, extracellular brown, nucleus teal, Golgi olive, vacuole purple).
- Intensity ramp: `0.15 + 0.85 × (expr / P95)^0.60`.
- Reactions without a GPR or with no expression: `#e7ebef` grey.

1045 of 1622 reactions are painted (all with a numeric expression score).

## Files & scripts

```
rna_seq_integration/
├── src/
│   ├── stage5_pathway_analysis.py     (classifier + enrichment + reporter + condition matrix)
│   └── build_expression_map.py        (expression-painted Escher HTML)
├── outputs/
│   ├── pathway_assignment.tsv         (per-reaction pathway assignment; auditability)
│   ├── pathway_summary.tsv            (one row per pathway with expr + enrichment stats)
│   ├── pathway_enrichment.tsv         (Fisher OR + p + BH q per pathway)
│   ├── pathway_condition_matrix.tsv   (pathway × condition × O2 flux totals)
│   ├── reporter_metabolites.tsv       (Patil-Nielsen z scores)
│   └── biomass_pathway_corr.tsv       (pathway ↔ biomass Pearson r)
└── reports/
    ├── PATHWAY_ANALYSIS.md            (this document)
    └── PATHWAY_ANALYSIS.docx          (formatted for manuscripts)
```

All raw output TSVs are linked as pill-buttons directly under the Pathway
breakdown section of the RNA-seq tab.

## References

- Patil KR, Nielsen J. *Uncovering transcriptional regulation of metabolism
  by using metabolic network topology.* PNAS 2005; 102:2685-9.
- Blazier AS, Papin JA. *Integration of expression data in genome-scale
  metabolic network reconstructions.* Front Physiol 2012; 3:299.
- Machado D, Herrgård M. *Systematic evaluation of methods for integration
  of transcriptomic data into constraint-based models of metabolism.*
  PLoS Comput Biol 2014; 10:e1003580.
