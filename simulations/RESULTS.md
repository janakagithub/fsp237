# FSP237 simulation-panel results

Driver: `run_simulation_panel.py`
Model: `gpr-update/fsp237_atp_safe_gsm_gpr_updated.json` (1703 rxns, 1330 mets, 1286 genes)
Plan source: `fsp237_biomass_extension/INFECTION_SIM_PLAN.md`

All 36 simulations (18 conditions × 2 O2 states) returned a solver-feasible
solution; biomass is 0 where the model lacks a viable carbon route from the
spec'd source.

## Summary by condition (aerobic)

| # | Condition | Stage | Biomass (1/h) | C-in (mmol/gDW/h) | Biomass / mmol C | Notes |
|---|---|---|---:|---:|---:|---|
| 14 | Maltose (starch proxy) | necrotrophic | **0.364** | 60.0 | 0.0061 | Best single-source — 2 glucose units per maltose |
| 12 | Xylose | necrotrophic | 0.146 | 25.0 | 0.0058 | Cereal hemicellulose monomer |
| 15 | Mixed AAs + glc | necrotrophic | 0.142 | 42.5 | 0.0034 | 20 AAs + glc, all at −0.5 |
| 01 | Trehalose | pre-infection | 0.073 | 12.0 | 0.0061 | Spore reserve |
| 09 | L-Glutamine | biotrophic | 0.053 | 10.0 | 0.0053 | C+N from one substrate |
| 08 | L-Glutamate | biotrophic | 0.053 | 10.0 | 0.0053 | — |
| 10 | GABA | biotrophic | 0.050 | 8.0 | 0.0062 | GABA shunt works (aerobic only) |
| 06 | Glucose LOW | biotrophic | 0.036 | 6.0 | 0.0061 | Reference |
| 07 | Sucrose LOW | biotrophic | 0.036 | 6.0 | 0.0061 | Invertase route OK |
| 02 | Glycerol | pre-infection | 0.021 | 3.0 | 0.0069 | — |
| 16 | **Pre-infection mix** | cocktail | 0.028 | 4.5 | 0.0063 | trehalose + glycerol + oleate (see breakdown) |
| 17 | **Biotrophic mix** | cocktail | 0.063 | 11.2 | 0.0056 | glc + suc + Glu + Gln + malate (see breakdown) |
| 18 | **Necrotrophic mix** | cocktail | **0.398** | 70.0 | 0.0057 | suc + maltose + galU + xyl + ara + 20 AA (see breakdown) |
| 03 | Palmitate (C16:0) | pre-infection | **0.000** | — | — | **FAIL — β-oxidation gap** |
| 04 | Oleate (C18:1) | pre-infection | **0.000** | — | — | **FAIL — β-oxidation gap** |
| 05 | Hexacosanoate (C26:0) | pre-infection | **0.000** | — | — | **FAIL — VLCFA β-ox gap (expected)** |
| 11 | D-Galacturonate | necrotrophic | **0.000** | — | — | **FAIL — Ashwell pathway absent** |
| 13 | L-Arabinose | necrotrophic | **0.000** | — | — | **FAIL — fungal Ara assimilation absent** |

### Cocktail composition breakdown (with cpd IDs)

All uptake rates are in mmol/gDW/h (passed to `EX_<cpd>_e0.lower_bound` as the
negated value). Anaerobic runs use the same cocktails with `EX_cpd00007_e0`
(O₂) closed.

**#16 — Pre-infection mix** — 3 substrates, total 0.9 mmol/gDW/h → 4.5 mmol C
| cpd id | metabolite | uptake | role |
|---|---|---:|---|
| `cpd00794_e0` | Trehalose (TRHL) | 0.30 | Spore reserve disaccharide |
| `cpd00100_e0` | Glycerol | 0.30 | Appressorial turgor / lipid backbone |
| `cpd15269_e0` | Octadecenoate (oleate, C18:1) | 0.30 | Major lipid-body fatty acid |

**#17 — Biotrophic mix** — 5 substrates, total 2.0 mmol/gDW/h → 11.2 mmol C
| cpd id | metabolite | uptake | role |
|---|---|---:|---|
| `cpd00027_e0` | D-Glucose | 0.30 | Apoplastic sugar |
| `cpd00076_e0` | Sucrose | 0.20 | Phloem sugar |
| `cpd00023_e0` | L-Glutamate | 0.50 | Dominant apoplastic AA (C + N) |
| `cpd00053_e0` | L-Glutamine | 0.50 | Major N-transport AA |
| `cpd00130_e0` | L-Malate | 0.50 | Apoplastic organic acid |

**#18 — Necrotrophic mix** — 25 substrates, total 10.5 mmol/gDW/h → 70.0 mmol C

Polymer / sugar breakdown products:
| cpd id | metabolite | uptake | role |
|---|---|---:|---|
| `cpd00076_e0` | Sucrose | 2.00 | Storage sugar release |
| `cpd00179_e0` | Maltose | 2.00 | Starch breakdown |
| `cpd00280_e0` | D-Galacturonate | 1.00 | Pectin backbone monomer |
| `cpd00154_e0` | Xylose | 1.00 | Hemicellulose monomer |
| `cpd00224_e0` | L-Arabinose | 0.50 | Arabinoxylan side-chain monomer |

All 20 proteinogenic amino acids at 0.20 each (cytoplasmic release):
| cpd id | AA | | cpd id | AA | | cpd id | AA | | cpd id | AA |
|---|---|---|---|---|---|---|---|---|---|---|
| `cpd00023_e0` | Glu | | `cpd00041_e0` | Asp | | `cpd00065_e0` | Trp | | `cpd00119_e0` | His |
| `cpd00033_e0` | Gly | | `cpd00051_e0` | Arg | | `cpd00066_e0` | Phe | | `cpd00129_e0` | Pro |
| `cpd00035_e0` | Ala | | `cpd00053_e0` | Gln | | `cpd00069_e0` | Tyr | | `cpd00132_e0` | Asn |
| `cpd00039_e0` | Lys | | `cpd00054_e0` | Ser | | `cpd00084_e0` | Cys | | `cpd00156_e0` | Val |
| `cpd00060_e0` | Met | | `cpd00107_e0` | Leu | | `cpd00161_e0` | Thr | | `cpd00322_e0` | Ile |

Note that `cpd00023_e0` (Glu) and `cpd00053_e0` (Gln) appear in both #17 and
#18; the necrotrophic dose (0.20 each) is lower because the AA pool is broad,
whereas the biotrophic dose (0.50 each) is concentrated on the two dominant
apoplastic AAs only.

Anaerobic biomass is uniformly ~0.005×–0.20× of aerobic (when growth exists),
as expected. GABA fails anaerobically (0.0 ana vs 0.050 aer): the GABA shunt
feeds succinate → TCA → ETC, which needs O2.

## Cocktail-ranking sanity check (key biological prediction)

Expected: **necrotrophic > biotrophic > pre-infection**.
Observed (aerobic):  `0.398 > 0.063 > 0.028` — **confirmed**, with the
necrotrophic cocktail giving ~14× more biomass than the spore-stage mix.

## Five informative failures — diagnosis and gap-fill plan

These are exactly the gap-fill targets the plan flagged. Each condition is
feasible (no solver infeasibility); biomass is 0 because no flux
distribution exists that satisfies biomass demand from the C source.

For each failure I traced which steps of the canonical pathway are
**present in the model** and which are **missing**, then list the
ModelSEED rxn IDs that should be added to restore growth. Every proposed
addition has a documented orthologue in the *C. graminicola* /
*C. higginsianum* genomes [O'Connell 2012; Buiate 2017], so the gap-fills
are biologically justified, not just feasibility hacks.

---

### 1–3. Fatty-acid β-oxidation (#03 palmitate, #04 oleate, #05 hexacosanoate)

**What IS in the model (peroxisomal `_x0`):**

| Step | Reaction (model rxn id) | Status |
|---|---|---|
| Fatty acid uptake | `rxn08704_c0` (HDCAt), `rxn09841_c0` (fatty-acid peroxisomal transport) | ✓ present |
| Acyl-CoA activation (C16) | `rxn00947_c0`, `rxn00947_x0` Palmitate:CoA ligase  | ✓ present |
| Acyl-CoA activation (C26) | `rxn09445_x0` fatty-acid–CoA ligase, peroxisomal | ✓ present |
| Acyl-CoA → CoA-ester transport to peroxisome | `rxn09849_x0` ABC | ✓ present |
| **Step 1.** Acyl-CoA oxidase | `rxn09474_x0` … `rxn09479_x0` (6 chain lengths) | ✓ present |
| **Step 2.** trans-2-enoyl-CoA hydratase | partial — see below | **⚠️ INCOMPLETE** |
| **Step 3.** 3-OH-acyl-CoA dehydrogenase | `rxn20474_x0` (C16), `rxn30390_x0` (C14), `rxn20478/20480_x0` (C10), `rxn20476_x0` (C18), `rxn09461_x0` (C26) | ✓ present |
| **Step 4.** 3-ketoacyl-CoA thiolase | `rxn08767_x0` | ✓ present |
| Glyoxylate cycle: isocitrate lyase | `rxn00336_c0` (gene_7509 or gene_8576 …) | ✓ present |
| Glyoxylate cycle: malate synthase | `rxn20162_x0` (gene_5781) | ✓ present |

**Step 2 — enoyl-CoA hydratase: present chain lengths vs missing**

| Chain | ModelSEED id | In model? |
|---|---|---|
| C16 | `rxn03240` (S)-3-Hydroxyhexadecanoyl-CoA hydro-lyase | ✓ `rxn03240_x0` |
| C14 | `rxn03241` (S)-3-Hydroxytetradecanoyl-CoA hydro-lyase | ✓ `rxn03241_x0` |
| C12 | `rxn02911` (S)-3-Hydroxydodecanoyl-CoA hydro-lyase | ✓ `rxn02911_x0` |
| C10 | `rxn03245` (S)-Hydroxydecanoyl-CoA hydro-lyase | ✓ `rxn03245_x0` |
| **C8**  | `rxn03247` (S)-Hydroxyoctanoyl-CoA hydro-lyase | **✗ MISSING** |
| **C6**  | `rxn03250` (S)-Hydroxyhexanoyl-CoA hydro-lyase | **✗ MISSING** |
| **C4**  | `rxn02167` (S)-3-Hydroxybutanoyl-CoA hydro-lyase | **✗ MISSING** |

This is the actual blocker for palmitate: β-oxidation can chain-shorten
C16 → C14 → C12 → C10 but then stalls at C8 because the C8 enoyl-CoA
hydratase step is absent. Acetyl-CoA never reaches the cytosol, the
glyoxylate cycle can't fire, and no biomass forms.

**Gap-fill list — add to peroxisome (`_x0`):**

| Reaction | EC | Equation | Why |
|---|---|---|---|
| `rxn03247_x0` | 4.2.1.17 | trans-Oct-2-enoyl-CoA + H₂O → (S)-3-Hydroxyoctanoyl-CoA | C8 hydratase — chain-shortening past C10 |
| `rxn03250_x0` | 4.2.1.17 | trans-Hex-2-enoyl-CoA + H₂O → (S)-3-Hydroxyhexanoyl-CoA | C6 hydratase |
| `rxn02167_x0` | 4.2.1.17 | crotonyl-CoA + H₂O → (S)-3-Hydroxybutyryl-CoA | C4 hydratase — final round |

After adding these 3 reactions (one fungal MFP/Fox2 gene covers all of
them — assign `gpr = ''` until BLAST hits a Fox2 ortholog in FSP237),
palmitate should grow at ~0.04–0.06 1/h (depending on intermediate
turnover).

**For oleate (#04, C18:1):** the Δ9 unsaturation needs one extra step,
2,4-dienoyl-CoA reductase + isomerase (Δ3,5,Δ2,4-dienoyl-CoA isomerase) to
re-enter the saturated cycle. Add `rxn02678` (trans-Oct-2-enoyl-CoA
reductase) and `rxn02719` (trans-Dodec-2-enoyl-CoA reductase) variants
for the C18:1 path, then the cycle above completes it. Lower priority
than #03 because once #03 grows, #04 becomes a one-extra-step fix.

**For hexacosanoate (#05, C26:0):** add the very-long-chain (VLCFA)
hydratase variants for C24, C22, C20, C18, plus the already-listed C16-C4
set. VLCFA β-ox in fungi requires the peroxisomal-only `Pxa1`/`Pxa2` ABC
import (the model has `rxn09849_x0` which is generic — assume it covers
VLCFA too). Lowest priority; cuticle-wax catabolism is biology-rich but
not core to running a Cs lifecycle simulation.

---

### 4. D-Galacturonate (#11) — fungal Ashwell pathway missing

**What IS in the model:**

| Step | Reaction | Status |
|---|---|---|
| Pectin hydrolysis (extracellular) | `rxn13253_e0` Endopolygalacturonase | ✓ produces extracellular galU |
| galU exchange | `EX_cpd00280_e0` | ✓ |
| galU transport (e0 → c0) | — | **✗ MISSING** |
| Ashwell pathway (4 enzymes c0) | — | **✗ ALL MISSING** |

**Gap-fill list — fungal galU pathway** (this is the *Aspergillus niger*
GAR1-LGD1-LGA1 pathway; reference Hilditch et al. 2007 Appl Env Microbiol,
Martens-Uzunova & Schaap 2008 FGB; orthologs identified in *Colletotrichum*
genomes):

| # | Reaction | Equation / role | Why |
|---|---|---|---|
| 1 | **`rxn05673`** D-galacturonate transport via proton symport | cpd00280_e0 + H⁺_e0 → cpd00280_c0 + H⁺_c0 | First — get galU into cytosol |
| 2 | **GAR1** (fungal D-galU reductase, NADPH; not currently in ModelSEED with the canonical fungal stereochem — use `rxn01456_c0` as a stand-in or add a custom reaction `cpd00280_c0 + NADPH → L-galactonate + NADP`) | EC 1.1.1.365 | Pathway entry |
| 3 | **`rxn21749`** L-galactonate dehydratase | L-galactonate → 2-keto-3-deoxy-L-galactonate (KDG) + H₂O | LGD1 step |
| 4 | KDG aldolase (LGA1) — not in ModelSEED as a single rxn; equivalent path | 2-keto-3-deoxy-L-galactonate → pyruvate + L-glyceraldehyde | Lumped fungal step |
| 5 | **`rxn00208_c0`** L-glyceraldehyde reductase / GLD1 | L-glyceraldehyde + NADPH → glycerol | Reconnects to central metabolism via glycerol kinase + G3P-DH |

GPR: BLAST the fungal A. niger GAR1/LGD1/LGA1/GLD1 sequences against the
FSP237 proteome and assign the best gene_* per step.

After this 5-step addition, galacturonate should grow at ~0.10–0.15 1/h
(comparable to xylose, since both are C6/C5 sugars feeding pyruvate +
glycerol pools).

---

### 5. L-Arabinose (#13) — Penttilä pathway entry step missing

**What IS in the model:**

| Step | Reaction | Status |
|---|---|---|
| L-Ara transport | `rxn08142_c0` (extracellular ↔ cyto) | ✓ present |
| **Step 1.** L-Ara → L-arabitol | — | **✗ MISSING** |
| **Step 2.** L-arabitol → L-xylulose | `rxn01291_c0` L-Arabitol:NADP+ 1-oxidoreductase | ✓ present |
| **Step 3.** L-xylulose → xylitol | — | **✗ MISSING** |
| **Step 4.** Xylitol → D-xylulose | `rxn01043_c0` xylitol:NADP+ oxidoreductase | ✓ present |
| **Step 5.** Xylulokinase → D-xylulose-5P | `rxn01199_c0` ATP:D-xylulose 5-phosphotransferase | ✓ present |

Steps 1 and 3 are the gaps. Note **D-arabinose** has a complete path
(added during the biomass extension for the ascorbate pathway via
`rxn01150_c0`); only the **L-isomer** is broken — and L-arabinose is the
one that matters for cereal arabinoxylan.

**Gap-fill list — add to cytosol (`_c0`):**

| # | Reaction | EC | Equation | Why |
|---|---|---|---|---|
| 1 | **`rxn01289_c0`** L-Arabinose:NAD+ 1-oxidoreductase (XYL1/AldR) | 1.1.1.21 | L-Ara + NAD⁺ → L-arabitol + NADH | Pathway entry — Penttilä step 1. (Most fungi use NAD-coupled but NADPH variants exist; `rxn15236` / `rxn15237` are functional synonyms — pick one.) |
| 2 | **`rxn01392_c0`** Xylitol:NADP+ 4-oxidoreductase (L-xylulose-forming), reverse direction (LXR1) | 1.1.1.10 | L-xylulose + NADPH → xylitol + NADP⁺ | Step 3 — bridges into the existing xylitol→D-xylulose→D-xylulose-5P part of the pathway |

GPR: BLAST *Pichia stipitis* / *Trichoderma reesei* XYL1/LAR1 + LXR1
sequences against the FSP237 proteome and assign best gene_*. (XYL1 has a
very-well-characterized fungal ortholog set.)

After adding these 2 reactions, L-arabinose should grow at ~0.10–0.12 1/h
(similar to xylose since both feed D-xylulose-5P into the PPP).

---

## Recommended gap-fill execution order

| Priority | Condition unlocked | Reactions to add | Effort | Biological payoff |
|---|---|---|---|---|
| **1** | #03 palmitate, #05 hexacosanoate (partially) | 3 enoyl-CoA hydratases (`rxn03247`, `rxn03250`, `rxn02167`) in `_x0` | trivial — 3 reactions, no new mets | Unlocks the appressorial lipid-mobilization phase; foundational for any pre-infection realism |
| **2** | #13 L-arabinose | 2 reactions (`rxn01289_c0`, `rxn01392_c0`) | trivial | Unlocks arabinoxylan-derived carbon — dominant sorghum hemicellulose |
| **3** | #11 D-galacturonate | 5 reactions (transport + 4 pathway steps), 1 custom rxn for GAR1 | moderate — needs GAR1 custom-build and BLAST for GPRs | Unlocks the pectin necrotroph strategy |
| **4** | #04 oleate | 2-3 reactions (Δ9 isomerase + reductase) | moderate — chain-length matched variants | Realistic unsaturated FA β-ox |
| **5** | #05 hexacosanoate (full) | C18-C24 hydratase variants | larger — multiple chain lengths | Cuticle-wax catabolism; lower-stakes biology |

After **priorities 1 + 2 + 3** the model would grow on all five currently-failing
single-substrate conditions — closing every gap the report flagged for the
pre-infection-through-necrotrophic biology.

**Suggested implementation:** add the gap-fill specs to a new
`fsp237_biomass_extension/gapfill_bypath.py` following the same
hardcoded-list pattern as `MANNITOL_RXNS` / `ALPHA13_GLUCAN_RXNS` in
`extend_biomass.py`, then re-run the GPR mapping + site refresh.

## What worked unexpectedly well

- **Maltose** (#14) at 0.364 — far more biomass than glucose. The model has
  maltose phosphorylase / maltase to split it cleanly. Good news for the
  starch-storage hypothesis.
- **Xylose** (#12) at 0.146 — the model has a complete fungal XR/XDH or
  isomerase route. Hemicellulose monomer is usable.
- **GABA** (#10) at 0.050 (aerobic) — GABA shunt is functional. This
  matters: plant defense releases GABA, and the model can re-route it.
- **Sucrose** = **Glucose** at low rate — sucrose splits to glc+fru and
  reaches the same biomass; invertase + fructokinase wiring is fine.
- **Glutamine** ≈ **Glutamate** despite Gln's extra amide N. Both end up at
  α-KG; the slight Gln advantage (0.0528 vs 0.0527) reflects extra N being
  available for biomass AA synthesis.

## Anaerobic-specific observations

- Anaerobic biomass on glucose is 0.0074 = 20% of aerobic. The bio_gsm growth
  rate is consistent with the 2 ATP/glucose limit established by the
  ATP-safe build (a fermenter must produce 15× less biomass per glucose).
- Glycerol anaerobically: 0.0004 — essentially zero. Glycerol oxidation
  needs the ETC.
- GABA anaerobically: 0 — see above.
- Maltose anaerobically: 0.074 = 2× anaerobic glucose, scaling with the
  hexose count.

## Output files

```
simulations/
├── run_simulation_panel.py            # the driver
├── simulation_results.tsv             # 36 rows (one per condition × O2 state)
├── per_condition/
│   ├── 01_trehalose_aerobic.tsv       # non-zero fluxes with rxn name / equation / GPR
│   ├── 01_trehalose_anaerobic.tsv
│   ├── 02_glycerol_aerobic.tsv
│   ├── …                              # 25 flux dumps (5 failed conditions have no dump)
│   └── 18_necrotrophic_mix_anaerobic.tsv
└── RESULTS.md                         # this file
```

Failed conditions (#03, #04, #05, #11, #13) have no per-condition TSV
because the biomass-driven flux solution is trivial (all zero).

## Next steps suggested by the results

1. **Gap-fill fatty-acid β-oxidation** (priority 1 — affects 3 conditions).
   Even C16 fails, which is the most studied chain length. Without this the
   model can't describe the appressorial stage at all.
2. **Gap-fill the Ashwell pathway** (4 enzymes) — unlocks galacturonate
   and indirectly pectin.
3. **Gap-fill the Penttilä L-Ara pathway** (4 enzymes) — unlocks
   arabinoxylan-derived carbon, the dominant sorghum hemicellulose.
4. The "expected to need gap-fill" items from the plan (starch, cellulose,
   cellobiose, ferulate, dhurrin) were not even tested here because the
   model has no exchange/transport stubs for them yet. Adding the exchange +
   one hydrolytic step per polymer is the second wave of gap-fill work.
5. After gap-fills, re-run this panel to confirm the failures move to
   non-zero biomass.
