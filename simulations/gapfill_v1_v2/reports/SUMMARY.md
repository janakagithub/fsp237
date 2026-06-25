# FSP237 gap-fill — Versions 1 & 2 summary

**Built:** 2026-06-22
**Source model:** `gpr-update/fsp237_atp_safe_gsm_gpr_updated.json`
**Pipeline:** `simulations/gapfill_v1_v2/`
**Biomass coefficients:** unchanged (per user instruction — only reactions
and GPRs are touched; the biomass reaction is identical to the source).

---

## Versioning trail

| Version | File | What's in it |
|---|---|---|
| **source** | `gpr-update/fsp237_atp_safe_gsm_gpr_updated.json` | 1703 rxns / 1330 mets / 1286 genes; GPR-clean, biomass-extended |
| **V1** | `models/fsp237_gapfilled_Version1_noGenes.json` | source + **28 new reactions** (no GPRs) → 1731 rxns / 1347 mets / 1286 genes |
| **V2** | `models/fsp237_gapfilled_Version2_gapfill_genes_integrated.json` | V1 + GPRs on the 19 confidently-mapped reactions → 1731 rxns / 1347 mets / **1293 genes** |
| **V3** | `models/fsp237_gapfilled_Version3_dedup_noGenes.json` | V1 with exact-duplicate cleanup (33 dup groups → 34 rxns removed; `bio1` protected) → **1697** rxns / 1347 mets / 1286 genes |
| **V4** | `models/fsp237_gapfilled_Version4_dedup_genes_integrated.json` | V2 with the same dedup pass, GPRs merged across duplicates → **1697** rxns / 1347 mets / 1293 genes |
| **V5** | `models/fsp237_gapfilled_Version5_dirlock_noGenes.json` | V3 with β-oxidation + Ashwell direction-locked to degradation only (34 reactions tightened) → **1697** rxns / 1293 genes (unchanged from V3) |
| **V6** | `models/fsp237_gapfilled_Version6_dirlock_genes_integrated.json` | V4 with the same direction lock applied → **1697** rxns / 1293 genes (unchanged from V4) |
| **V9** | `models/fsp237_gapfilled_Version9_vlcfa_complete_noGenes.json` | V5 + (a) complete VLCFA chain C26→C24→C22→C20→C18 (13 new rxns incl. transport, biosynthesis-consistent cpd IDs), (b) dropped 88 reactions in g0/n0/v0 (zero-flux compartments) → **1621** rxns / 1268 mets / 1267 genes |
| **V10** | `models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json` | V6 + the same VLCFA + compartment cleanup → **1622** rxns / 1268 mets / 1274 genes (currently published on the site) |

Each version is preserved on disk; V3 was built from V1 (not from source),
V4 from V2, V5 from V3, V6 from V4 -- so the gap-fill / gene-integration /
dedup / direction-lock states are cleanly separable and reversible.
V3 ≡ V4 except for GPRs; V5 ≡ V6 except for GPRs.

## Exact-duplicate dedup (V3 / V4)

User identified `ALCD2ir_c0` ≡ `rxn00543_c0` (ethanol DH) as an example of
duplicates not caught by earlier passes. The `build_v3v4_dedup.py` script
scans the whole model for exact-stoichiometry (compartment-aware, direction-
normalized) duplicates and applies a precedence policy:

1. ModelSEED `rxn*` IDs beat BiGG-style aliases (e.g., keep `rxn00543_c0`,
   drop `ALCD2ir_c0`).
2. Among ModelSEED IDs, the lower numeric suffix wins (matches the
   `BuildMinimalFSP237Model.ipynb` dedup convention -- e.g. `rxn00102` beats
   `rxn30392`).
3. `rxn` prefix beats `frxn` (canonical beats fix).
4. Alphabetical tie-break.

**GPRs are not lost:** every duplicate's `gene_reaction_rule` is union-merged
into the kept reaction's GPR (deduplicating gene tokens, preserving order
of first appearance). Bounds are widened to the superset of all
duplicates' bounds.

Results:
- **33 duplicate groups** after `bio1` protection (was 34; `bio1` is
  removed from its group so both it and `rxn00062_c0` survive)
- **34 reactions removed** (one group has 3 members: `rxn09524_c0`,
  `frxn08975_c0`, `rxn09524_m0` -- the cross-compartment NADH-DH triple)
- **9 groups had a GPR transfer** (the kept canonical picked up extra
  genes from the dropped duplicates)
- **14 groups had widened bounds** (e.g., dropped duplicate was reversible
  while kept was forward-only, so the kept becomes reversible)
- **0 metabolites removed** (a duplicate's metabolites are also used by the
  kept canonical, so none became orphans)

**Protected reactions (not subject to dedup):**

| Rxn ID | Why protected |
|---|---|
| `bio1` | ATP-yield objective used by `extend_biomass.py`, `build_atp_safe_site.py`, `build_escher_maps.py`, and the simulation panel for the canonical 30/2 ATP/glucose test. Stoichiometrically identical to `rxn00062_c0` (NGAM) but kept as a separate role-marked reaction with no GPR. Without protection, the dedup pass would collapse `bio1 → rxn00062_c0` and break every ATP-yield helper that hardcodes the `bio1` objective. Per user instruction 2026-06-23. |

The `PROTECTED_IDS` set in `build_v3v4_dedup.py` is the place to add more
such role-marker reactions if other duplicates need preservation.

## Direction lock (V5 / V6)

User instruction (2026-06-23): the gap-filled β-oxidation cycle and the
Ashwell D-galacturonate pathway should be **unidirectional in the
degradation direction only**, so the optimizer cannot run them in reverse
to synthesize palmitoyl-CoA from acetyl-CoA, or galacturonate/L-galactonate
from glycerol -- those are not biologically realistic for a fungal
pathogen.

`build_v5v6_dirlock.py` reads V3/V4 and tightens 34 reactions' bounds
according to the per-rxn `DIRECTION_LOCKS` spec, where each entry says
which sign of flux in the ModelSEED-written equation corresponds to
degradation (some rxns are written synthesis-forward, some degradation-
forward -- the spec accounts for both). Locking happens by:

- **forward-locked**: `bounds = (0, upper)` -- only positive flux allowed
- **reverse-locked**: `bounds = (lower, 0)` -- only negative flux allowed

### What's locked

| Block | Reactions | Lock direction | Net effect |
|---|---|---|---|
| Beta-ox step 1 (acyl-CoA oxidase / FAD-DH) | 11 rxns (C26/C18/C16/C14/C12/C10/C8/C6/C4) | forward | acyl-CoA → enoyl-CoA only |
| Beta-ox step 2 (enoyl-CoA hydratase) | 9 rxns | mostly reverse | enoyl + H₂O → 3-OH only (ModelSEED writes them backwards for most chains) |
| Beta-ox step 3 (3-OH-acyl-CoA DH) | 9 rxns | mixed (3 forward + 6 reverse) | 3-OH → 3-oxo only |
| Beta-ox step 4 (chain-shortening thiolase) | 7 rxns | reverse | 3-oxo + CoA → AcCoA + (n-2) only |
| Ashwell D-galU pathway | 5 rxns | mixed | galU import-only; galU → L-galactonate → KDG → pyruvate + glycerol only |

### What's NOT locked (deliberately)

- **Peroxisomal cofactor shuttles (tx_atp_xc, tx_ppi_xc, tx_coa_xc,
  tx_nad_xc, tx_nadh_xc, tx_h_xc)** -- transport stubs that need to be
  reversible for cycle cofactor regeneration.
- **Acyl-CoA ligases (`rxn00947_c0/_x0`, `rxn09445_x0`, `rxn05736_x0`)** --
  fatty-acid activation isn't part of the chain-shortening cycle per se.
- **Penttilä L-arabinose pathway** -- user explicitly out of scope for the
  direction-lock pass; remains reversible (steps 2 + 3 are `<=>`).
- **Pectin → galU extracellular hydrolysis (`rxn13253_e0`)** -- already
  forward-only in source.

### Verification

- All 36 panel conditions × 2 O₂ states produce bit-identical biomass
  (within 1e-9) to V3/V4. Locks didn't open or close any productive path;
  they only blocked unused reverse paths.
- Spot-check on palmitate-only condition: the C10→C8 thiolase
  (`rxn02680_x0`) carries flux of -0.72 (reverse = degradation), all
  other locked rxns at zero or correct sign.
- Spot-check on necrotrophic cocktail: GAR1 = -1.0 (reverse = NADPH-
  reduction direction = degradation), LGD1/LGA1/GLD1 = +1.0 (forward =
  degradation). All Ashwell enzymes engaged correctly with no reverse
  flow.

Full per-rxn bound audit in `reports/v5_direction_locks.tsv`.

Notable dedup decisions (full audit in `reports/v3_dedup_log.tsv` and
`reports/v4_dedup_log.tsv`):

| Kept | Dropped | Reason |
|---|---|---|
| `rxn00543_c0` ethanol-DH | `ALCD2ir_c0` | ModelSEED beats BiGG alias |
| `rxn09524_c0` NADH-DH | `frxn08975_c0`, `rxn09524_m0` | `rxn` > `frxn`; cross-compartment duplicate consolidated |
| `rxn00102_c0` carbonic anhydrase | `rxn30392_c0` | Lower numeric suffix wins |
| `rxn00003_m0` PDC | `rxn33164_m0` | Lower numeric suffix wins |
| `rxn02680_x0` C10→C8 thiolase | `rxn19988_x0` | Same biology, lower numeric suffix; V1 added rxn02680, source had rxn19988 -- they collapse |
| `rxn02804_x0` C16→C14 thiolase | `rxn19979_x0` | Same biology; gene_434, gene_5400, gene_9875 transferred from rxn19979 |
| `rxn03243_x0` C12→C10 thiolase | `rxn19982_x0` | Same biology; gene_1935, gene_5400 transferred |
| `rxn06510_x0` C14→C12 thiolase | `rxn19986_x0` | Same biology; gene_1935, gene_5400 transferred |
| `bio1` ATP demand | `rxn00062_c0` | Lower-ID `bio1` was the existing source canonical |
| `rxn09882_c0` | `rxn09882_m0` | Cross-compartment near-identical; consolidated to cytosol |
| (... 24 more documented in the TSV log) | | |

**Note on the β-ox thiolases:** four of my V1 additions
(`rxn02804`, `rxn06510`, `rxn03243`, `rxn02680`) had functional duplicates
in the source model (`rxn19979`, `rxn19986`, `rxn19982`, `rxn19988`). The
dedup pass consolidated each pair into the lower-numeric-suffix variant,
keeping my V1 additions and merging in the GPRs (gene_1935, gene_5400,
gene_434, gene_9875) that the source map already had on those thiolases.
Net effect on V4: those 4 V1-added thiolases now have richer GPRs than
they did in V2.

**Validation:** V3 and V4 panel results are bit-identical to V1 and V2
(within solver noise) across all 18 conditions × 2 O2 states. No biomass
change, no ATP-yield change. Confirmed by `test_gapfilled_model.py v3`
and `v4` against the 36-row baseline.

---

## Phase 1 — gap-fill reactions added in V1 (no GPRs)

All stoichiometry pulled from ModelSEED; per-rxn details in
`reports/v1_added_reactions.tsv`.

### A. Peroxisomal β-oxidation closure — 22 reactions in `_x0`

The existing model had `_x0` acyl-CoA oxidases + 3-OH-DHs for C18/C16/C14/C12/C10
and C26, but was missing the C8/C6/C4 cycle steps and 6 of the 7 chain-shortening
thiolases. Added per user requirement to keep β-ox in the peroxisome (no
cytosolic lumping).

| Step | Chain | Reaction (rxn_id) | EC |
|---|---|---|---|
| Acyl-CoA dehydrogenase (oxidase) | C8 | `rxn02679_x0` | 1.3.3.6 |
| | C6 | `rxn03251_x0` | 1.3.3.6 |
| | C4 | `rxn00868_x0` | 1.3.99.2 |
| Enoyl-CoA hydratase | C8 | `rxn03247_x0` | 4.2.1.17 |
| | C6 | `rxn03250_x0` | 4.2.1.17 |
| | C4 | `rxn02167_x0` | 4.2.1.55 |
| 3-OH-acyl-CoA dehydrogenase | C8 | `rxn03246_x0` | 1.1.1.35 |
| | C6 | `rxn03249_x0` | 1.1.1.35 |
| | C4 (acetoacetyl-CoA red.) | `rxn03861_x0` | 1.1.1.36 |
| 3-Ketoacyl-CoA thiolase (chain-shortening) | C14→C12 | `rxn06510_x0` | 2.3.1.16 |
| | C12→C10 | `rxn03243_x0` | 2.3.1.16 |
| | C10→C8 | `rxn02680_x0` | 2.3.1.16 |
| | C8→C6 | `rxn03248_x0` | 2.3.1.16 |
| | C6→C4 | `rxn00874_x0` | 2.3.1.16 |
| | C4 → 2 acetyl-CoA | `rxn00178_x0` | 2.3.1.9 |
| Peroxisomal cofactor shuttles | ATP/AMP antiport | `tx_atp_xc` | — |
| | PPi export | `tx_ppi_xc` | — |
| | CoA equilibration | `tx_coa_xc` | — |
| | NAD equilibration | `tx_nad_xc` | — |
| | NADH equilibration | `tx_nadh_xc` | — |
| | H+ equilibration | `tx_h_xc` | — |

(C16→C14 thiolase `rxn02804_x0` was already present in the source model
under a different name pattern — skipped on add as `already_present`.)

The shuttle reactions are required for flux: without ATP/CoA/NAD transport
into and out of the peroxisome, the cycle's cofactor pools are closed loops
and the model would still show zero β-ox flux even with all chain-length
enzymes in place. Diagnosed by max-producing each cofactor in palmitate
conditions and observing 0 / 0 / 0 yields pre-shuttles.

### B. Penttilä L-arabinose pathway closure — 2 reactions in `_c0`

The model already had L-Ara → L-arabitol (`rxn01291_c0`) and the downstream
xylitol → D-xylulose → D-xylulose-5P part. Missing the middle:

| # | Reaction | Substrate → product | EC |
|---|---|---|---|
| 1 | `rxn01391_c0` L-arabitol dehydrogenase (NAD) | L-arabitol → L-xylulose | 1.1.1.12 |
| 2 | `rxn33066_c0` L-xylulose reductase (NADPH) | L-xylulose → xylitol | 1.1.1.10 |

### C. Fungal D-galacturonate (Ashwell) pathway — 5 reactions in `_c0`

Pectin can already be hydrolyzed to galU extracellularly (`rxn13253_e0`),
but the cytosolic Aspergillus-type Ashwell pathway was absent. Added:

| # | Reaction | Function | EC |
|---|---|---|---|
| 1 | `rxn05673_c0` galU/H+ symporter | e0 → c0 transport | 2.A.1.14.- |
| 2 | `rxn07491_c0` D-galU NADPH-reductase (GAR1) | galU → L-galactonate | 1.1.1.365 |
| 3 | `rxn21749_c0` L-galactonate dehydratase (LGD1) | L-galactonate → 2-keto-3-deoxy-L-galactonate | 4.2.1.146 |
| 4 | `rxn21750_c0` KDG aldolase (LGA1) | KDG → pyruvate + L-glyceraldehyde | 4.1.2.55 |
| 5 | `rxn09954_c0` L-glyceraldehyde reductase (GLD1) | L-glyceraldehyde + NADPH → glycerol | 1.1.1.21 |

---

## Phase 1 — V1 test results vs source

Simulation panel results in `reports/v1_simulation_results.tsv`. Aerobic
biomass for the conditions that were previously blocking:

| Condition | Source biomass | V1 biomass | Change |
|---|---:|---:|---|
| #03 Palmitate (C16) | 0.000 | **0.114** | ✓ now grows; full peroxisomal β-ox closed |
| #04 Oleate (C18:1) | 0.000 | **0.128** | ✓ now grows; β-ox cycle + Δ3-Δ2 isom path |
| #05 Hexacosanoate (C26) | 0.000 | 0.000 | ✗ still fails — VLCFA needs more chain variants |
| #11 D-Galacturonate | 0.000 | **0.149** | ✓ now grows; full Ashwell |
| #13 L-Arabinose | 0.000 | **0.140** | ✓ now grows; Penttilä closed |
| All other 31 (condition × O2) | unchanged | unchanged | (within 1e-4) |

**Cocktails also improved:**
- #16 Pre-infection mix: 0.028 → **0.069** (oleate now contributes)
- #17 Biotrophic mix: 0.063 → 0.063 (unchanged; was C-limited not gap-limited)
- #18 Necrotrophic mix: 0.398 → **0.446** (galU + L-Ara now contribute)

Solver status: every condition `optimal`; ATP yield 30/2 preserved.

---

## Phase 2 + 3 — Colletotrichum gene candidates → FSP237 BLAST

Pipeline in `find_candidates.py`. 47 candidate queries pulled from
`gpr-update/C_higgensium.gbff` (RefSeq C. higginsianum IMI 349063,
CH63R_ locus tags) by scoring each CDS against the gap-fill rxn's EC set
+ product-name keywords; top 1-2 candidates retained per rxn; full BLAST
against the existing FSP237 BLAST DB at `gpr-update/blast_db/fsp237`.

Quality thresholds (matching the original GPR-overhaul pipeline):
**pident ≥ 30, qcov ≥ 50, e-value ≤ 1e-10**. Confidence bands:
- **high**: pident ≥ 60 AND qcov ≥ 70
- **medium**: pident ≥ 45 AND qcov ≥ 60
- **low**: passes thresholds but below medium
- **none**: no candidate / no passing hit

Full per-rxn detail: `blast/v1_rxn_to_gene_mapping.tsv` (final picks),
`blast/v1_blast_results.tsv` (raw BLAST), `candidates/v1_query_provenance.tsv`
(which CH63R protein was the query, with product description).

**Summary of mapping outcomes** (28 V1 gap-fill reactions):

| Confidence | Count | Reactions |
|---|---|---|
| high | 22 | All 17 β-ox chain steps + 5 of 6 shuttles + L-Ara step 1 + Ashwell GAR1 + GLD1 |
| medium | 2 | `rxn01391_c0` Penttilä step 2 (gene_8784, pid 48%); `rxn05673_c0` galU transport (gene_5858, pid 46%) |
| low | 0 | — |
| **none** | **4** | `rxn33066_c0` L-xylulose reductase; `rxn07491_c0` GAR1; `rxn21749_c0` LGD1; `rxn21750_c0` LGA1 |

For the 4 unmapped reactions: the C. higginsianum genome's RefSeq /product=
annotations don't carry the canonical Aspergillus-style enzyme names
("LXR1", "LGD1", "LGA1"), so the keyword filter missed them. The genes
almost certainly exist in C. higginsianum (these pathways are pan-fungal
in pectinolytic Ascomycetes per Buiate 2017 BMC Genomics), but they need
a different identification strategy: BLAST FSP237 directly against the
Aspergillus niger or Trichoderma reesei reference protein for each enzyme,
or against the InterPro signature.

---

## Phase 4 — Version 2 with integrated GPRs

`build_v2_integrate_genes.py` reads `v1_rxn_to_gene_mapping.tsv` and writes
the picked `fsp237_gene` onto each V1 reaction as a single-gene GPR.

**Conservative posture:**
- High + medium confidence: GPR assigned.
- Low or none: GPR left blank, recorded in audit trail.
- **Shuttle reactions**: GPR left blank by default (`FLAG_SHUTTLES_UNASSIGNED = True`
  in script) because the BLAST keyword filter used generic terms (e.g.,
  "MFS transporter") that don't pinpoint the actual peroxisomal membrane
  carrier. Flip the flag if you want the BLAST picks integrated.

**Audit trail:** `reports/v2_gpr_assignments.tsv` (one row per V1 gap-fill
rxn with: decision, gpr assigned, source CH63R locus tag, source product
description, source strain, BLAST metrics, reason).

**Assignment counts:** 19 reactions assigned a GPR; 10 left unassigned
(4 Ashwell/Penttilä without C. higg hit + 6 shuttle stubs flagged for review).

V2 model invariants:
- 1731 reactions (same as V1)
- 1293 genes (1286 + 7 newly referenced: gene_2621, gene_2953, gene_8837,
  gene_11751, gene_11942, gene_4118, gene_9875). The genes for shuttles
  + medium-confidence picks would add a few more if we relaxed the
  conservative posture.
- Biomass aer/ana = 0.1821 / 0.0370 (unchanged from source)
- ATP/glc aer/ana = 30 / 2 (unchanged)
- All Panel results identical to V1 (genes don't affect FBA flux, only
  enzyme-knockout simulations).

---

## Reactions still without confident gene assignment

These 10 are worth a second pass with curated Aspergillus/Trichoderma
reference sequences (the genus-level pathway literature is good):

| Rxn | Why no assignment | Recommended next step |
|---|---|---|
| `rxn33066_c0` L-xylulose reductase | No C. higg /product= keyword match | BLAST FSP237 with *T. reesei* LXR1 (UniProt G0RBI5) |
| `rxn07491_c0` GAR1 (D-galU reductase) | Only one weak C. higg "Gar1" hit -- that's a different protein (RNA-binding) | BLAST FSP237 with *A. niger* GaaA (CAK39597.1) |
| `rxn21749_c0` LGD1 (L-galactonate dehydratase) | No keyword match | BLAST FSP237 with *A. niger* GaaB (CAK40253.1) |
| `rxn21750_c0` LGA1 (KDG aldolase) | No keyword match | BLAST FSP237 with *A. niger* GaaC (CAK46208.1) |
| `tx_atp_xc` Ant1-like peroxisomal ATP/AMP | Generic "carrier" hit, low specificity | BLAST FSP237 with *S. cerevisiae* Ant1p (YPR128C) |
| `tx_ppi_xc` PPi export | No specific transporter | Likely doesn't need a gene (passive transport) |
| `tx_coa_xc` CoA equilibration | No specific transporter | Likely no single gene; leave unassigned |
| `tx_nad_xc` NAD equilibration | No specific transporter | Leave unassigned (likely passive) |
| `tx_nadh_xc` NADH equilibration | No specific transporter | Leave unassigned (likely passive) |
| `tx_h_xc` H+ equilibration | No specific transporter | Leave unassigned (likely passive) |

---

## Deliverables (all under `simulations/gapfill_v1_v2/`)

```
build_v1_gapfill.py                              # Phase 1 builder
build_v2_integrate_genes.py                      # Phase 4 builder
find_candidates.py                               # Phase 2+3 candidate finder + BLAST
test_gapfilled_model.py <model> <suffix>         # Re-run simulation panel against any model

models/
  fsp237_gapfilled_Version1_noGenes.json
  fsp237_gapfilled_Version2_gapfill_genes_integrated.json

candidates/
  v1_query_candidates.faa                        # 47 C. higginsianum protein queries
  v1_query_provenance.tsv                        # per-query metadata

blast/
  v1_blast_results.tsv                           # raw blastp output (183 lines)
  v1_rxn_to_gene_mapping.tsv                     # final picks per rxn

reports/
  v1_added_reactions.tsv                         # 28 V1 reactions with status, EC, equation, reason
  v1_simulation_results.tsv                      # 36 simulations on V1
  v1_per_condition/                              # per-condition flux dumps for V1
  v2_simulation_results.tsv                      # 36 simulations on V2 (identical to V1; GPRs don't change FBA)
  v2_per_condition/                              # per-condition flux dumps for V2
  v2_gpr_assignments.tsv                         # audit trail of every gap-fill rxn -> GPR decision
  SUMMARY.md                                     # this file
```

## How to revert

Both V1 and V2 are full standalone JSON models. To go back to the V1 model
(gap-fill but no genes), load `models/fsp237_gapfilled_Version1_noGenes.json`.
To go back to the source (pre-gap-fill), load
`../../gpr-update/fsp237_atp_safe_gsm_gpr_updated.json`. Nothing in the
gap-fill pipeline modifies upstream files.
