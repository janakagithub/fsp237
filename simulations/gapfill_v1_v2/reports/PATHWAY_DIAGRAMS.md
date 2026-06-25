# Degradation pathways gap-filled in the FSP237 model

**Model versions:** all pathways below are present and fully wired in
**V1, V2, V3, V4** (`simulations/gapfill_v1_v2/models/`).
**V3/V4** are the recommended deliverable (post-dedup).

For each pathway:
- a step diagram from substrate to a central-metabolism handoff
- a table of every reaction in the diagram, with its ID, **human-readable equation**, EC number, and whether it was **(existing)** in the source model or **(V1 added)** in the gap-fill pass
- the GPRs (V2/V4 only) where assigned

In equations, compartment suffixes are dropped for readability except where
crossing membranes (`_e0` extracellular, `_c0` cytosol, `_x0` peroxisome).
"`L-Lyxitol`" and "`L-Lyxulose`" are the ModelSEED canonical names for the
KEGG `L-arabitol` and `L-xylulose` -- same molecules, different naming
conventions.

---

## Pathway 1 — Peroxisomal β-oxidation of saturated long-chain fatty acids

**Goal:** any peroxisomal acyl-CoA (C26, C18, C16, C14, C12, C10, C8, C6, C4)
is degraded one round of `oxidase → hydratase → 3-OH-DH → thiolase`,
shortening by 2 carbons per round, producing 1 acetyl-CoA per round, until
two acetyl-CoAs remain. Acetyl-CoA exits to the glyoxylate cycle / TCA via
the existing peroxisomal `rxn20162_x0` malate synthase and the cytosolic
`rxn00336_c0` isocitrate lyase.

### Membrane transport + activation (entry to the cycle)

Palmitate is the example below; identical paths exist for the other chain
lengths via `rxn09445_x0` (C26 ligase) and the `rxn08448/08449/08451/08457_x0`
chain-shortening ligases.

```
   Palmitate (extracellular)
   cpd00214_e0
        │     HDCAt_c0   (facilitated diffusion, irreversible)
        ▼
   Palmitate (cytosol)
   cpd00214_c0
        │     rxn08704_c0   (Palmitate transport c0 → x0)
        ▼
   Palmitate (peroxisome)
   cpd00214_x0
        │     rxn00947_x0   (Palmitate:CoA ligase; ATP + CoA → AMP + PPi + Palmitoyl-CoA)
        ▼
   Palmitoyl-CoA (peroxisome)
   cpd00134_x0  ← cytosolic Palmitoyl-CoA also enters peroxisome
                   via rxn09849_x0 (fatty-acyl-CoA ABC transporter)
        │
        ▼
   [Enter the β-oxidation cycle below]
```

### β-oxidation cycle (one round per chain length)

```
                       ╭───────────────────────────────────────────╮
                       │                                           │
   Acyl-CoA (Cn)       │                                           │
   ─────────►  STEP 1: acyl-CoA oxidase / dehydrogenase            │
   + O2 or FAD         │   produces H2O2 (oxidase) or FADH2 (DH)   │
                       ▼                                           │
   trans-2-Enoyl-CoA (Cn)                                          │
                       │                                           │
                       ▼  STEP 2: enoyl-CoA hydratase  (+ H2O)     │
   (S)-3-Hydroxyacyl-CoA (Cn)                                      │
                       │                                           │
                       ▼  STEP 3: 3-OH-acyl-CoA dehydrogenase      │
   + NAD                  produces NADH                            │
                       ▼                                           │
   3-Oxoacyl-CoA (Cn)                                              │
                       │                                           │
                       ▼  STEP 4: 3-ketoacyl-CoA thiolase          │
   + CoA                  produces 1 acetyl-CoA per round          │
                       ▼                                           │
   Acyl-CoA (Cn-2)  +  Acetyl-CoA  ────────────► out to glyoxylate │
                       │                            cycle / TCA    │
                       │                                           │
                       ╰─── repeat with new (Cn-2) substrate ──────╯
```

### Chain-length-resolved enzyme table

The model now has all four cycle steps for every chain length from C16
down to C4 (C18 + C26 also covered). Steps marked **(V1 added)** were
inserted by the gap-fill pass; **(existing)** were already in the source
GPR-updated model. **(absorbed ⊕ rxnXXXXX)** means a V1-added reaction
was merged with a source-model duplicate during the V3/V4 dedup pass,
keeping the V1 ID and absorbing the source rxn's GPRs.

| Chain | Step | Reaction ID | Human-readable equation | EC | Source |
|---|---|---|---|---|---|
| **C26** | 1. oxidase | `rxn09474_x0` | O₂ + Hexacosanoyl-CoA → H₂O₂ + trans-Hexacos-2-enoyl-CoA | 1.3.3.6 | existing |
|  | 2. hydratase | `rxn09473_x0` | H₂O + trans-Hexacos-2-enoyl-CoA ⇌ (S)-3-Hydroxyhexacosyl-CoA | 4.2.1.17 | existing |
|  | 3. 3-OH-DH | `rxn09461_x0` | NAD + (S)-3-Hydroxyhexacosyl-CoA ⇌ NADH + H⁺ + 3-Oxohexacosyl-CoA | 1.1.1.35 | existing |
|  | 4. thiolase to C24 | — | — | 2.3.1.16 | **gap (not yet filled)** |
| **C18** | 1. oxidase | `rxn09475_x0` | O₂ + Stearoyl-CoA → H₂O₂ + (2E)-Octadecenoyl-CoA | 1.3.3.6 | existing |
|  | 1. acyl-CoA DH (FAD) | `rxn08053_x0` | FAD + Stearoyl-CoA ⇌ FADH₂ + (2E)-Octadecenoyl-CoA | 1.3.99.- | existing |
|  | 2. hydratase | `rxn08413_x0` | (S)-3-Hydroxyoctadecanoyl-CoA ⇌ H₂O + (2E)-Octadecenoyl-CoA | 4.2.1.17 | existing |
|  | 3. 3-OH-DH | `rxn20476_x0` | NADH + H⁺ + 3-Oxostearoyl-CoA ⇌ NAD + (S)-3-Hydroxyoctadecanoyl-CoA | 1.1.1.35 | existing |
|  | 4. thiolase C18→C16 | `rxn08767_x0` | CoA + 3-Oxostearoyl-CoA → Acetyl-CoA + Palmitoyl-CoA | 2.3.1.16 | existing |
| **C16** | 1. oxidase | `rxn09476_x0` | O₂ + Palmitoyl-CoA → H₂O₂ + (2E)-Hexadecenoyl-CoA | 1.3.3.6 | existing |
|  | 2. hydratase | `rxn03240_x0` | (S)-3-Hydroxyhexadecanoyl-CoA ⇌ H₂O + (2E)-Hexadecenoyl-CoA | 4.2.1.17 | existing |
|  | 3. 3-OH-DH | `rxn20474_x0` | NADH + H⁺ + 3-Oxopalmitoyl-CoA ⇌ NAD + (S)-3-Hydroxyhexadecanoyl-CoA | 1.1.1.35 | existing |
|  | 4. thiolase C16→C14 | `rxn02804_x0` | Acetyl-CoA + Myristoyl-CoA ⇌ CoA + 3-Oxopalmitoyl-CoA | 2.3.1.16 | V1 added (⊕ rxn19979) |
| **C14** | 1. oxidase | `rxn09477_x0` | O₂ + Myristoyl-CoA → H₂O₂ + (2E)-Tetradecenoyl-CoA | 1.3.3.6 | existing |
|  | 2. hydratase | `rxn03241_x0` | (S)-3-Hydroxytetradecanoyl-CoA ⇌ H₂O + (2E)-Tetradecenoyl-CoA | 4.2.1.17 | existing |
|  | 3. 3-OH-DH | `rxn30390_x0` | NADH + H⁺ + 3-Oxotetradecanoyl-CoA ⇌ NAD + (S)-3-Hydroxytetradecanoyl-CoA | 1.1.1.35 | existing |
|  | 4. thiolase C14→C12 | `rxn06510_x0` | Acetyl-CoA + Lauroyl-CoA ⇌ CoA + 3-Oxotetradecanoyl-CoA | 2.3.1.16 | V1 added (⊕ rxn19986) |
| **C12** | 1. oxidase | `rxn09478_x0` | O₂ + Lauroyl-CoA → H₂O₂ + (2E)-Dodecenoyl-CoA | 1.3.3.6 | existing |
|  | 1. acyl-CoA DH (FAD) | `rxn02720_x0` | FAD + Lauroyl-CoA ⇌ FADH₂ + (2E)-Dodecenoyl-CoA | 1.3.99.- | existing |
|  | 2. hydratase | `rxn02911_x0` | (S)-3-Hydroxydodecanoyl-CoA ⇌ H₂O + (2E)-Dodecenoyl-CoA | 4.2.1.17 | existing |
|  | 3. 3-OH-DH | `rxn20478_x0` | NADH + H⁺ + 3-Oxododecanoyl-CoA ⇌ NAD + (S)-3-Hydroxydodecanoyl-CoA | 1.1.1.35 | existing |
|  | 4. thiolase C12→C10 | `rxn03243_x0` | Acetyl-CoA + Decanoyl-CoA ⇌ CoA + 3-Oxododecanoyl-CoA | 2.3.1.16 | V1 added (⊕ rxn19982) |
| **C10** | 1. oxidase | `rxn09479_x0` | O₂ + Decanoyl-CoA → H₂O₂ + (2E)-Decenoyl-CoA | 1.3.3.6 | existing |
|  | 2. hydratase | `rxn03245_x0` | (S)-Hydroxydecanoyl-CoA ⇌ H₂O + (2E)-Decenoyl-CoA | 4.2.1.17 | existing |
|  | 3. 3-OH-DH | `rxn20480_x0` | NADH + H⁺ + 3-Oxodecanoyl-CoA ⇌ NAD + (S)-Hydroxydecanoyl-CoA | 1.1.1.35 | existing |
|  | 4. thiolase C10→C8 | `rxn02680_x0` | Acetyl-CoA + Octanoyl-CoA ⇌ CoA + 3-Oxodecanoyl-CoA | 2.3.1.16 | V1 added (⊕ rxn19988) |
| **C8** | 1. acyl-CoA DH (FAD) | `rxn02679_x0` | FAD + Octanoyl-CoA ⇌ FADH₂ + (2E)-Octenoyl-CoA | 1.3.3.6 | **V1 added** |
|  | 2. hydratase | `rxn03247_x0` | (S)-Hydroxyoctanoyl-CoA ⇌ H₂O + (2E)-Octenoyl-CoA | 4.2.1.17 | **V1 added** |
|  | 3. 3-OH-DH | `rxn03246_x0` | NAD + (S)-Hydroxyoctanoyl-CoA ⇌ NADH + H⁺ + 3-Oxooctanoyl-CoA | 1.1.1.35 | **V1 added** |
|  | 4. thiolase C8→C6 | `rxn03248_x0` | Acetyl-CoA + Hexanoyl-CoA ⇌ CoA + 3-Oxooctanoyl-CoA | 2.3.1.16 | **V1 added** |
| **C6** | 1. acyl-CoA DH (FAD) | `rxn03251_x0` | FAD + Hexanoyl-CoA ⇌ FADH₂ + (2E)-Hexenoyl-CoA | 1.3.3.6 | **V1 added** |
|  | 2. hydratase | `rxn03250_x0` | (S)-Hydroxyhexanoyl-CoA ⇌ H₂O + (2E)-Hexenoyl-CoA | 4.2.1.17 | **V1 added** |
|  | 3. 3-OH-DH | `rxn03249_x0` | NAD + (S)-Hydroxyhexanoyl-CoA ⇌ NADH + H⁺ + 3-Oxohexanoyl-CoA | 1.1.1.35 | **V1 added** |
|  | 4. thiolase C6→C4 | `rxn00874_x0` | Acetyl-CoA + Butyryl-CoA ⇌ CoA + 3-Oxohexanoyl-CoA | 2.3.1.16 | **V1 added** |
| **C4** | 1. acyl-CoA DH (NAD) | `rxn00868_x0` | NAD + Butyryl-CoA ⇌ NADH + H⁺ + Crotonyl-CoA | 1.3.99.2 | **V1 added** |
|  | 2. hydratase | `rxn02167_x0` | (S)-3-Hydroxybutanoyl-CoA ⇌ H₂O + Crotonyl-CoA | 4.2.1.55 | **V1 added** |
|  | 3. 3-OH-DH (NADP) | `rxn03861_x0` | NADP + (S)-3-Hydroxybutyryl-CoA ⇌ NADPH + H⁺ + Acetoacetyl-CoA | 1.1.1.36 | **V1 added** |
|  | 4. terminal thiolase | `rxn00178_x0` | 2 Acetyl-CoA ⇌ CoA + Acetoacetyl-CoA | 2.3.1.9 | **V1 added** |

### Peroxisomal cofactor shuttles (all V1 added)

Without these the cofactor pools in the peroxisome are closed loops and
no β-ox flux can carry. Each is a reversible passive-diffusion stub
(refine to specific antiporters in a future curation step).

| Reaction ID | Equation | Role |
|---|---|---|
| `tx_atp_xc` | ATP_c0 + AMP_x0 ⇌ ATP_x0 + AMP_c0 | ATP/AMP antiport (Ant1-like) — funds acyl-CoA activation |
| `tx_ppi_xc` | PPi_x0 ⇌ PPi_c0 | Export PPi produced by acyl-CoA ligase |
| `tx_coa_xc` | CoA_c0 ⇌ CoA_x0 | Free-CoA equilibration |
| `tx_nad_xc` | NAD_c0 ⇌ NAD_x0 | NAD import for 3-OH-DH steps |
| `tx_nadh_xc` | NADH_c0 ⇌ NADH_x0 | NADH export (via cytosolic redox shuttle) |
| `tx_h_xc` | H⁺_c0 ⇌ H⁺_x0 | Proton equilibration |

### Acetyl-CoA exit to central metabolism (existing)

```
   Acetyl-CoA (peroxisome)  ─────► malate synthase  (rxn20162_x0)
   cpd00022_x0                     Acetyl-CoA + Glyoxalate → Malate + CoA
                                       │
                                       ▼
                                   Malate (x0 / c0)
                                       │
                                       │  cytosolic isocitrate lyase
                                       │  (rxn00336_c0):
                                       │    Isocitrate → Succinate + Glyoxalate
                                       │
                                       └─► glyoxylate cycle → gluconeogenesis → biomass
```

---

## Pathway 2 — Penttilä L-arabinose assimilation

**Goal:** extracellular L-arabinose is reduced to L-arabitol, oxidized to
L-xylulose, reduced again to xylitol, then oxidized to D-xylulose and
phosphorylated to D-xylulose-5-phosphate, which feeds the pentose
phosphate pathway. This is the canonical fungal pathway
(Penttilä et al.; *Pichia stipitis* / *Trichoderma reesei*); it is
distinct from the bacterial L-arabinose isomerase route.

```
   L-Arabinose (extracellular)
   cpd00224_e0
        │     rxn08142_c0   (L-Ara transport e0 ⇌ c0)
        ▼
   L-Arabinose (cytosol)
   cpd00224_c0
        │     rxn01291_c0  (L-Ara reductase, NADPH)        [existing — was Step 1]
        ▼
   L-Arabitol (L-Lyxitol in ModelSEED naming)
   cpd00417_c0
        │     rxn01391_c0  (L-arabitol DH, NAD)            [V1 added — Step 2]
        ▼
   L-Xylulose (L-Lyxulose in ModelSEED naming)
   cpd00261_c0
        │     rxn33066_c0  (L-xylulose reductase, NADPH)   [V1 added — Step 3]
        ▼
   Xylitol
   cpd00306_c0
        │     rxn01043_c0  (Xylitol DH, NADP)              [existing — Step 4]
        ▼
   D-Xylulose
   cpd00154_c0
        │     rxn01199_c0  (Xylulokinase, ATP)             [existing — Step 5]
        ▼
   D-Xylulose-5-phosphate  ──► Pentose phosphate pathway → biomass
```

### Reaction table

| # | Reaction ID | Equation | EC | Source |
|---|---|---|---|---|
| 0 | `rxn08142_c0` | L-Arabinose_e0 ⇌ L-Arabinose_c0 | (transport) | existing |
| 1 | `rxn01291_c0` | NADPH + H⁺ + L-Arabinose_c0 → NADP + L-Arabitol_c0 | 1.1.1.21 | existing |
| 2 | `rxn01391_c0` | NAD + L-Arabitol_c0 ⇌ NADH + H⁺ + L-Xylulose_c0 | 1.1.1.12 | **V1 added** |
| 3 | `rxn33066_c0` | NADP + Xylitol_c0 ⇌ NADPH + H⁺ + L-Xylulose_c0 | 1.1.1.10 | **V1 added** |
| 4 | `rxn01043_c0` | NADPH + H⁺ + D-Xylose_c0 → NADP + Xylitol_c0 (runs reverse for Penttilä) | 1.1.1.9 | existing |
| 5 | `rxn01199_c0` | ATP + D-Xylulose_c0 → ADP + H⁺ + D-Xylulose-5-phosphate | 2.7.1.17 | existing |

**Diagnostic test:** L-arabinose as sole carbon source — biomass went
from 0 (source model) to **0.140 mmol/gDW/h aerobic** in V1+ after
adding steps 2 and 3.

---

## Pathway 3 — Fungal D-galacturonate (Ashwell) degradation

**Goal:** pectin (already hydrolyzed extracellularly to galacturonate
by the existing `rxn13253_e0`) is transported into the cytosol, reduced
to L-galactonate, dehydrated to 2-keto-3-deoxy-L-galactonate (KDG), cleaved
to pyruvate + L-glyceraldehyde, and the glyceraldehyde reduced to glycerol
which enters central metabolism via the existing glycerol kinase. This is
the *Aspergillus niger* / *Trichoderma reesei* GAR1-LGD1-LGA1-GLD1 pathway
(Hilditch et al. 2007; Martens-Uzunova & Schaap 2008).

```
   Pectin (extracellular)                    [existing extracellular hydrolysis]
   cpd11601_e0
        │     rxn13253_e0   (Endopolygalacturonase)
        │     H2O + Pectin → D-Galacturonate
        ▼
   D-Galacturonate (extracellular)
   cpd00280_e0
        │     rxn05673_c0   (galU/H⁺ symporter)              [V1 added]
        │     H⁺_e0 + D-Galacturonate_e0 ⇌ H⁺_c0 + D-Galacturonate_c0
        ▼
   D-Galacturonate (cytosol)
   cpd00280_c0
        │     rxn07491_c0   (D-galU reductase / GAR1)         [V1 added]
        │     NADPH + H⁺ + D-galU ⇌ NADP + L-Galactonate
        ▼
   L-Galactonate
   cpd14659_c0
        │     rxn21749_c0   (L-galactonate dehydratase / LGD1) [V1 added]
        │     L-Galactonate → H₂O + 2-keto-3-deoxy-L-galactonate
        ▼
   2-keto-3-deoxy-L-galactonate (KDG)
   cpd23364_c0
        │     rxn21750_c0   (KDG aldolase / LGA1)             [V1 added]
        │     KDG ⇌ Pyruvate + L-Glyceraldehyde
        ▼
   Pyruvate + L-Glyceraldehyde
   cpd00020_c0   cpd01605_c0
                   │     rxn09954_c0   (L-glyceraldehyde reductase / GLD1) [V1 added]
                   │     NADPH + H⁺ + L-Glyceraldehyde ⇌ NADP + Glycerol
                   ▼
                Glycerol (cpd00100_c0)
                   │                     [existing glycerol kinase + G3P-DH]
                   ▼
                Glycerol-3P → Dihydroxyacetone-P → glycolysis
   Pyruvate ─────────────────────────────► TCA / fermentation
```

### Reaction table

| # | Reaction ID | Equation | EC | Source |
|---|---|---|---|---|
| 0 | `rxn13253_e0` | H₂O + Pectin_e0 → D-Galacturonate_e0 | 3.2.1.15 | existing |
| 1 | `rxn05673_c0` | H⁺_e0 + D-Galacturonate_e0 ⇌ H⁺_c0 + D-Galacturonate_c0 | (transport) | **V1 added** |
| 2 | `rxn07491_c0` | NADPH + H⁺ + D-Galacturonate_c0 ⇌ NADP + L-Galactonate | 1.1.1.365 (GAR1) | **V1 added** |
| 3 | `rxn21749_c0` | L-Galactonate → H₂O + 2-keto-3-deoxy-L-galactonate | 4.2.1.146 (LGD1) | **V1 added** |
| 4 | `rxn21750_c0` | 2-keto-3-deoxy-L-galactonate ⇌ Pyruvate + L-Glyceraldehyde | 4.1.2.55 (LGA1) | **V1 added** |
| 5 | `rxn09954_c0` | NADPH + H⁺ + L-Glyceraldehyde ⇌ NADP + Glycerol | 1.1.1.21 (GLD1) | **V1 added** |

**Diagnostic test:** D-galacturonate as sole carbon source — biomass
went from 0 (source model) to **0.149 mmol/gDW/h aerobic** in V1+.

---

## Pathway summary

| Pathway | Reactions added in V1 | Was substrate growable before? | Now (V1/V3) |
|---|---:|---|---|
| Peroxisomal β-oxidation closure (C8/C6/C4 cycle + cofactor shuttles) | 16 + 6 shuttles | Palmitate 0.000, oleate 0.000 | Palmitate 0.114, oleate 0.128 |
| Penttilä L-arabinose | 2 | L-Ara 0.000 | L-Ara 0.140 |
| Fungal Ashwell galacturonate | 5 | galU 0.000 | galU 0.149 |
| **Total gap-fill reactions in V1/V3** | **29** | — | — |

## Gene-integrated versions (V2 / V4)

The V2 model adds single-gene GPRs from C. higginsianum BLAST hits to
16 of these new reactions (high-confidence hits). The V4 model adds the
dedup pass on top of V2, which merged in additional FSP237 gene tokens
on the chain-shortening thiolases (gene_1935, gene_5400, gene_434,
gene_9875, gene_4395) from absorbed `rxn199xx` duplicates. Per
`reports/v2_gpr_assignments.tsv` and `reports/v4_dedup_log.tsv`.
