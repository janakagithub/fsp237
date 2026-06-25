# FSP237 simulation plan: nutrient/carbon conditions across the *Colletotrichum sublineola* infection cycle on sorghum

**Model:** `fsp237_biomass_extension/fsp237_atp_safe_gsm_extended.json`
**Compiled:** 2026-06-22

---

## Section 1 — Summary

*Colletotrichum sublineola* (Cs) is a hemibiotrophic ascomycete that causes
anthracnose on *Sorghum bicolor*. Its infection cycle has three metabolically
distinct phases — surface/appressorial (no host access, lipid- and
trehalose-fueled), biotrophic primary hyphae in living cells (apoplastic
sugars and amino acids, low and competitive), and necrotrophic secondary
hyphae after cell lysis (host cytoplasm + cell-wall polysaccharides). Each
phase exposes the fungus to a different substrate panel, and a useful FBA
test regime needs to walk through all three, not just glucose.

The proposed plan is **15 single-source tests + 3 stage-mimicking cocktails**:

- **Pre-infection (4 conditions)**: trehalose, palmitate, oleate, hexacosanoate — exogenous-storage and cuticle-wax mobilization.
- **Biotrophic (5 conditions)**: glucose-low, sucrose-low, glutamate, GABA, glutamine — apoplastic substrates accessible during living-cell parasitism.
- **Necrotrophic (6 conditions)**: starch-derived maltose, sucrose, pectin → galacturonate, cellulose → cellobiose, hemicellulose → xylose, mixed amino acids — broad polymer breakdown.
- **3 cocktails** combining each stage's signature substrates.

Out of these, **5 conditions are expected to need gapfilling** (hexacosanoate β-oxidation, cellulose hydrolysis, hemicellulose/arabinoxylan branch sugars, starch hydrolysis to maltose, sorghum-specific dhurrin pathway). Each gap, if confirmed, is a *biologically real* metabolic capacity the model lacks — fixing it improves predictive power, not just feasibility.

All conditions use the existing extended model's `EX_cpd*_e0` exchanges; cpd IDs are given so each can be tested immediately by toggling the lower bound.

---

## Section 2 — Recommended nutrient/carbon conditions by infection stage

### A. Pre-infection / surface stage (0–24 h post-spore-deposit)

**Biological setting.** The spore lands on the waxy sorghum leaf surface,
germinates, and forms a melanized appressorium. The cuticle is hydrophobic
and effectively N-, P-, and labile-C-starved. Fungal nutrition is
**autonomous** in this phase: trehalose, mannitol, and lipid bodies in the
spore are mobilized, β-oxidation supplies acetyl-CoA, glyoxylate cycle
recycles C2 units into TCA intermediates, and glycerol accumulates to
generate appressorial turgor [Talbot 2003 ARM; Wang & Valent 2017].
Sorghum-leaf surface waxes are dominated by very-long-chain alkanes
(C29 and C31 predominate), primary alcohols (C26–C30), and free fatty
acids C24–C28 [Atkin & Hamilton 1982 J Nat Prod; Yates et al. 1991].
Cutin monomers (ω-hydroxy and dihydroxy C16/C18 acids) become available
once cutinase activity starts.

**Substrates the model should be able to use:**

| Substrate | Model cpd | Need gapfill? | Notes |
|---|---|---|---|
| Trehalose | `cpd00794_e0` | No | Endogenous reserve; biomass already includes intracellular TRHL |
| Glycerol | `cpd00100_e0` | No | Appressorial osmolyte / turgor |
| Palmitate (C16:0) | `cpd00214_e0` | No | Storage lipid hydrolysis product |
| Stearate (C18:0, `ocdca`) | `cpd01080_e0` | No | Storage lipid |
| Octadecenoate (C18:1, oleate) | `cpd15269_e0` | No | Major lipid-body FA |
| Octadecadienoate (C18:2, linoleate) | `cpd16888_e0` | No | Major lipid-body FA |
| Hexacosanoate (C26:0) | `cpd15240_e0` | **Possibly** | Cuticle wax FA; β-ox of very-long-chain may be incomplete |
| Cutin monomers (10,16-dihydroxypalmitate; 9,10,18-trihydroxystearate) | — | **Yes** | Not in model — would require cutinase output + ω-hydroxy-FA degradation |
| n-Alkanes (C25–C31) | — | **Yes** | Not in model — needs P450 alkane hydroxylase + ω-oxidation pathway |

### B. Infection / early biotrophic stage (24–72 h)

**Biological setting.** After penetration peg breaches the cuticle and cell
wall (via cutinases + minimal cell-wall enzymes — the biotrophic strategy
relies on *not* destroying the host yet [O'Connell et al. 2012 Nat Genet]),
an infection vesicle forms inside the first epidermal cell. Primary hyphae
(BIH) grow biotrophically through living cells, separated from the host
cytoplasm by an intact plasma membrane. Nutrients are limited to whatever
the host transporters secrete into the apoplast plus what the fungus can
extract via plant SWEET/STP sugar transporters that get manipulated by
effectors [Chen et al. 2010 Nature; Bezrutczyk et al. 2018 New Phytol].
The apoplast is **sugar-poor** (sub-mM glucose/fructose/sucrose), with low
levels of amino acids (Glu/Gln dominant in cereals), some organic acids
(malate, citrate), and the defense-amino-acid GABA which can be re-routed
as N/C by the pathogen [Solomon & Oliver 2002 Planta; Bolton 2009 MPMI].

**Stage-specific transcriptomic signal (from Cgr/Chig — closest evidence):**
during BIH stage, plant cell-wall degrading enzyme (PCWDE) genes are
**down-regulated** vs necrotrophic stage; secreted effector + secondary
metabolism genes are sharply up [O'Connell et al. 2012; Vargas et al. 2012
Plant Cell].

**Substrates the model should be able to use:**

| Substrate | Model cpd | Need gapfill? | Notes |
|---|---|---|---|
| Glucose (low, ≤1 mmol/gDW/h) | `cpd00027_e0` | No | Apoplastic; rate-limit not glucose-presence |
| Sucrose (low) | `cpd00076_e0` | No | Phloem sugar; needs invertase — already in iMM904 backbone |
| Fructose (low) | `cpd00082_e0` | No | Sucrose hydrolysis product |
| L-Glutamate | `cpd00023_e0` | No | Dominant cereal apoplast AA |
| L-Glutamine | `cpd00053_e0` | No | Dominant N transport AA |
| L-Asparagine | `cpd00132_e0` | No | Cereal N transport AA |
| L-Aspartate | `cpd00041_e0` | No | TCA-linked |
| GABA | `cpd00281_e0` | No | Host defense AA — fungal repurposing |
| L-Malate | `cpd00130_e0` | No | Apoplastic organic acid |
| Citrate | `cpd00137_e0` | No | Apoplastic |
| Methionine, Lys, Pro | `cpd00060`, `cpd00039`, `cpd00129_e0` | No | Minor apoplastic AAs |
| Dhurrin / p-hydroxybenzaldehyde (cyanogenic glucoside; sorghum-specific) | — | **Yes** | Sorghum-distinguishing defense compound; aglycone is toxic — Cs detoxification pathway is biologically relevant but completely absent from the model |

### C. Post-infection / necrotrophic stage (72 h+)

**Biological setting.** Cs switches from BIH to thicker, faster-growing
necrotrophic secondary hyphae (NSH). PCWDEs are massively up-regulated
[O'Connell 2012; Vargas 2012]. Host cells are killed and lyse, releasing
cytoplasm (free amino acids, sugars, organic acids, starch) and exposing
the cell wall. **Grass (Poaceae) cell walls** are Type II, with cellulose
microfibrils embedded in glucuronoarabinoxylan (GAX, dominant
hemicellulose), mixed-linkage β-(1,3;1,4)-glucan, low pectin (~5%), and
substantial ferulate cross-links between arabinose decorations on GAX
[Vogel 2008 Curr Opin Plant Biol; Carpita 1996 Annu Rev PMB]. Sorghum
also accumulates **starch** (esp. in stem and grain tissues — major C
reserve), free **proline** (osmotic adjustment in drought), and during
senescence/decay accumulates free **sucrose and hexoses** from polymer
breakdown.

**Substrates the model should be able to use:**

| Substrate | Model cpd | Need gapfill? | Notes |
|---|---|---|---|
| Sucrose (high) | `cpd00076_e0` | No | Major phloem/storage sugar |
| Maltose | `cpd00179_e0` | No (consumption); **Yes** for *production from starch* | Starch breakdown product |
| Starch (amylose/amylopectin) | — | **Yes** | Not in model as exchange; need α-amylase + glucoamylase to maltose/glucose |
| Cellobiose | — | **Yes** | Cellulose breakdown intermediate; need β-glucosidase + transport |
| Cellulose | — | **Yes** | Not feasible directly; expected to need cellobiohydrolase + transport stub |
| Xylose | `cpd00154_e0` | No (uptake); **Yes** for hemicellulose chain | Hemicellulose monomer |
| L-Arabinose / D-Arabinose | `cpd00224_e0` / `cpd00185_e0` | No | Arabinoxylan side-chain monomer; D-arabinose path was added in extension |
| D-Galacturonate | `cpd00280_e0` | No | Pectin backbone monomer |
| Pectin (homogalacturonan) | `cpd11601_e0` | No (exchange exists) | Should already produce galacturonate — verify in flux |
| 1,3-β-Glucan | `cpd11791_e0` | No | Cell wall + own wall material; verify uptake mode |
| Glycerol | `cpd00100_e0` | No | Lipid backbone from membrane lysis |
| Mixed amino acids (proteinogenic) | various `cpd000*_e0` | No | All 20 present as exchanges |
| Lecithin / phosphatidylcholine | `cpd11624_e0` | No | Membrane lipid release |
| Ferulate / p-coumarate (lignin-bound) | — | **Yes** | Sorghum cell wall is rich in hydroxycinnamates; degradation generates aromatic acids — *not in model* |

---

## Section 3 — Recommended simulation table

| # | Stage | Carbon source | Model cpd | Medium concept | Biological rationale | Expected behavior | Gapfill? | Notes |
|---|---|---|---|---|---|---|---|---|
| 1 | Pre-infection | Trehalose | `cpd00794_e0` | minimal + trehalose, –1 | Spore reserve mobilization | Growth at modest rate (~0.3–0.4× glc) | No | Tests trehalase + hexokinase route |
| 2 | Pre-infection | Glycerol | `cpd00100_e0` | minimal + glycerol, –1 | Appressorial osmolyte; also lipid backbone | Growth (gluconeogenesis + glycerol-3-P) | No | Confirms G3P dehydrogenase wired |
| 3 | Pre-infection | Palmitate (C16:0) | `cpd00214_e0` | minimal + palmitate, –1; glyoxylate ON | Lipid-body β-oxidation + glyoxylate cycle | Growth (requires ICL/MLS) | No | Diagnostic for glyoxylate shunt |
| 4 | Pre-infection | Oleate (C18:1) | `cpd15269_e0` | minimal + oleate, –1 | Major lipid-body FA | Growth; tests Δ9 desaturase reverse-flux | No | — |
| 5 | Pre-infection | Hexacosanoate (C26:0) | `cpd15240_e0` | minimal + C26 FA, –1 | Cuticle-wax FA | **Test:** may stall if very-long-chain β-ox missing | Possibly | Sorghum surface wax breakdown |
| 6 | Biotrophic | Glucose (LOW, –1) | `cpd00027_e0` | minimal, glc=−1 | Apoplastic limitation | Linear growth, low | No | Reference low-C condition |
| 7 | Biotrophic | Sucrose (LOW, –0.5) | `cpd00076_e0` | minimal, sucrose=−0.5 | Phloem sugar | Growth at ~2× glc-low (2 hexoses/sucrose) | No | Tests invertase + sugar mobilization |
| 8 | Biotrophic | L-Glutamate | `cpd00023_e0` | minimal + Glu –2 (C+N) | Apoplastic AA serving both | Growth via α-KG → TCA; N from amide | No | Should also satisfy N demand |
| 9 | Biotrophic | L-Glutamine | `cpd00053_e0` | minimal + Gln –2 (C+N) | Major N-transport AA | Growth, supplies 2N + α-KG | No | — |
| 10 | Biotrophic | GABA | `cpd00281_e0` | minimal + GABA –2 | Plant-defense AA repurposed | Growth via GABA shunt → succinate | No | If FAILS, model misses GABA shunt — important fix |
| 11 | Necrotrophic | Galacturonate (pectin) | `cpd00280_e0` | minimal + galU –5 | Pectin breakdown | Growth via Ashwell pathway | No | Reference necrotroph polymer C |
| 12 | Necrotrophic | Xylose | `cpd00154_e0` | minimal + Xyl –5 | Hemicellulose monomer | Growth via XR/XDH or isomerase route | No | Tests xylose assimilation |
| 13 | Necrotrophic | L-Arabinose | `cpd00224_e0` | minimal + Ara –5 | Arabinoxylan side-chain monomer | Growth via arabinose isomerase route | No | Cereal-specific |
| 14 | Necrotrophic | Maltose (starch proxy) | `cpd00179_e0` | minimal + maltose –5 | Starch breakdown product | Growth via maltase | No | Acceptable starch proxy without amylase |
| 15 | Necrotrophic | Mixed amino acids | 10–20 AA exchanges | rich-AA + glc=−0.5 | Cytoplasmic release | High biomass; N+C from AA | No | Cell-lysis simulation |
| 16 | Cocktail | **Pre-infection mix** | trehalose+glycerol+oleate at –0.3 each | starvation but lipid-fueled | Spore-stage realism | Growth at low rate | No | Should respire heavily |
| 17 | Cocktail | **Biotrophic mix** | glc=−0.3, suc=−0.2, glu=−0.5, gln=−0.5, mal=−0.5 | apoplast cocktail | Living-cell parasitism | Growth at intermediate rate | No | Tests AA + sugar integration |
| 18 | Cocktail | **Necrotrophic mix** | suc=−2, maltose=−2, galU=−1, xyl=−1, ara=−0.5, mixed AA=−0.2 each | host cytoplasm + cell wall | Lysed-cell realism | **Highest** biomass | No | Should outperform any single C |

---

## Section 4 — Highest-priority conditions to test first

Ranked by what each result tells you about the model and the biology:

1. **#6 Glucose-low (biotrophic reference)** — baseline. If this doesn't match the existing baseline (biomass ≈ 0.18 aerobic), something drifted.
2. **#3 Palmitate** — diagnostic for the glyoxylate cycle wiring. Cs MUST grow on lipids alone during the appressorial stage; if it can't, the model needs ICL/MLS curation.
3. **#1 Trehalose** — spore-reserve test. Should grow at roughly 2× hexose rate.
4. **#11 Galacturonate** — necrotrophic polymer breakdown; Ashwell pathway is a standard Colletotrichum capacity. Failure flags a CAZyme-downstream gap.
5. **#10 GABA** — biotrophic re-routing of host defense AA. Functional GABA shunt is a hallmark of phytopathogen success in living tissue.
6. **#17 Biotrophic cocktail** — first integrative test. Should reveal whether the model handles mixed limiting substrates correctly.
7. **#18 Necrotrophic cocktail** — should produce the highest biomass of any condition.
8. **#16 Pre-infection cocktail** — should produce the lowest biomass (correctly).
9. **#5 Hexacosanoate (C26)** — *failure here is expected and informative.* Tells you exactly which very-long-chain β-ox steps to add.
10. **#12, #13 Xylose / Arabinose** — verify hemicellulose-monomer assimilation works.

Run **1–8 first**; their outcomes will shape what gapfilling is worth doing for the lower-priority conditions.

---

## Section 5 — Key biological questions the simulations should answer

**Q1. Can the model support the metabolic transition from lipid-only (appressorial) to sugar-rich (necrotrophic)?**
This is the most fundamental hemibiotroph test. If conditions 1–4 grow at non-zero rate and 18 grows at the highest rate, the model can describe the full lifecycle in principle. If lipid-only conditions fail, the glyoxylate cycle wiring or β-oxidation is incomplete.

**Q2. Which host-derived substrates yield the highest biomass per C-atom?**
Sort conditions 6–15 by `biomass / C-uptake-rate`. Expected winners: sucrose, mixed-AA, glutamine (cofactor-of-2 because it's C+N). Expected laggards: galacturonate, arabinose. A surprise in this ranking points at a flux bottleneck.

**Q3. Does GABA grow Cs?**
Plant defense releases GABA; *successful* fungal pathogens reroute it via the GABA shunt (GABA → succinic semialdehyde → succinate). If condition #10 fails, the model is missing a documented Colletotrichum capacity — gapfill GABA transaminase + SSA dehydrogenase.

**Q4. Does the model expose the right starch-handling gap?**
Maltose alone (cond. #14) should work; full starch (no exchange) is expected to fail. Confirming this guides whether to add an α-amylase + glucoamylase chain or a starch-→-maltose lumped stub.

**Q5. Is there evidence the model needs cuticle-wax β-oxidation extension?**
Condition #5 (C26) tests this. Sorghum surface is rich in C26-C30 alkanes/FA. A robust Cs model should grow on C26 at non-zero rate (slow but positive). Failure → extend β-oxidation peroxisomal cycle by ≥4 rounds.

**Q6. Does the model show the expected biotrophic → necrotrophic biomass jump?**
Compare cocktails #17 vs #18. The literature describes a ~3–5× faster necrotrophic growth than biotrophic [O'Connell 2012]. If #18 is not noticeably higher than #17, mixed-substrate routing or NGAM may be miscalibrated.

**Q7. Where does the model break first?**
The 5 expected-gapfill conditions (#5, starch, cellulose, cellobiose, ferulate/dhurrin) are *intentionally* in the panel to expose biologically real gaps. Each failure is a curated `extend_biomass.py`-style follow-up.

---

## Evidence-quality footnote

- **Direct Cs–sorghum evidence**: Buiate et al. 2017 BMC Genomics (Cs vs Cgr comparative); Latunde-Dada & Lucas 2007 Mol Plant Pathol (Cs infection on sorghum); Baroncelli et al. 2014 Genome Announc 2:e00540-14 (Cs TX430BB genome). These confirm Cs is a hemibiotroph with expanded CAZyme repertoire.
- **Closest indirect** (other graminicolous *Colletotrichum*): Vargas et al. 2012 Plant Cell (*C. graminicola* on maize, full appressorial-to-necrotrophic biology); O'Connell et al. 2012 Nat Genet (*C. higginsianum* + *C. graminicola* in planta transcriptomes — defined the BIH/NSH transition framework).
- **Hemibiotrophy + appressorium metabolism** (genus-level / kingdom-level): Münch et al. 2008 Mycol Res; Talbot 2003 Annu Rev Microbiol; Wang & Valent 2017 Trends Microbiol.
- **Grass cell wall / sorghum surface chemistry**: Vogel 2008 Curr Opin Plant Biol; Carpita 1996 Annu Rev PMB; Atkin & Hamilton 1982 J Nat Prod; Yates et al. 1991 (sorghum wax).
- **Apoplast nutrient panel**: Solomon & Oliver 2002 Planta; Bolton 2009 MPMI; Bezrutczyk et al. 2018 New Phytol (SWEET transporters).
- **Sorghum-specific compounds**: dhurrin biosynthesis/distribution (Gleadow & Møller 2014 Annu Rev Plant Biol); arabinoxylan + ferulate cross-links (Burton & Fincher 2014 Curr Opin Plant Biol).

Where direct Cs evidence is sparse, the simulation plan defaults to the
*C. graminicola* + *C. higginsianum* paradigm (same hemibiotrophic
lifestyle, related genome scale, well-studied in planta transcriptomes).
This is the standard practice in the field and is the basis for the
sorghum-context choices above; any condition flagged "Yes" for gapfill
should ultimately be confirmed against Cs-specific in vitro carbon-source
growth data when available.
