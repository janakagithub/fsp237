# Proposed biomass edits — fsp237 / *Colletotrichum sublineola* FSP237

**Context.** Current `bio1` is `R_BIOMASS_SC5_notrace` from iMM904 (S. cerevisiae). FSP237 is a new, unpublished *C. sublineola* isolate; the closest published reference genome is **TX430BB** (Baroncelli et al. 2014, 46.75 Mb, 52.70% G+C, 12,699 genes). For exploration and v1 modeling, TX430BB is used as the proxy for FSP237.

**Coefficient units.** mmol per gDW. "Current" = value in iMM904 bio1. "Suggested" = literature-informed starting point; final values should be set after normalizing the total biomass to 1 g dry weight.

**Implementation priority.** Order: chitin → mannitol → ergosterol/zymosterol swap → α-1,3-glucan → Q9 swap → sphingolipid lump → nucleotide GC rebalance.

**Pathway column legend.**
- ✅ = in model already, no gap-fill needed
- 🆕 = new reaction(s) to add (ModelSEED IDs and equations given in the "MS rxn IDs + equations" column)
- ➖ = no biosynthesis change; coefficient-only

**ModelSEED equation conventions** in the "MS rxn IDs + equations" column: compartment suffixes (`_c0`, `_m0`, `_r0`) are omitted for brevity — every reaction is cytosolic unless noted. `cpd00067` = H⁺, `cpd00001` = H₂O, `cpd00009` = Pi, `cpd00012` = PPi.

---

## 1. Cell wall / structural polysaccharides

| Compound | Action | Metabolite cpd | Current | Suggested | Pathway | MS rxn IDs + equations | Evidence |
|---|---|---|---|---|---|---|---|
| **Chitin** | ADD to biomass | `cpd11683_c0` | 0 | **0.10–0.18** | ✅ complete | (no new rxn) — existing: `rxn09486_c0` GlcNAc-1-P + UTP ⇌ UDP-GlcNAc + PPi; `rxn09490_c0` cpd00175 → cpd00014 + cpd00067 + **cpd11683** | Bowman & Free 2006, *Bioessays* 28:799-808, PMID [16927300](https://pubmed.ncbi.nlm.nih.gov/16927300/), [DOI 10.1002/bies.20441](https://doi.org/10.1002/bies.20441). ssNMR on *A. fumigatus*/*A. nidulans* walls: Gautam et al. 2024, *Carb Polym*, PMID [39562136](https://pubmed.ncbi.nlm.nih.gov/39562136/) |
| **α-1,3-glucan** | ADD to biomass + gap-fill 1 reaction | `cpd12148_c0` (α-D-glucan polymer; already in ModelSEED, not in fsp237) | 0 | **0.05–0.20** (start ~0.08) | 🆕 1 rxn | **rxn15561** (EC 2.4.1.183) UDP-glucose:α-D-(1-3)-glucan 3-α-D-glucosyltransferase: cpd00026 (UDP-Glc) ⇌ cpd00014 (UDP) + **cpd12148** (α-1,3-glucan) | Fujikawa et al. 2012, *PLoS Pathog* 8:e1002882, PMID [22927818](https://pubmed.ncbi.nlm.nih.gov/22927818/), [DOI 10.1371/journal.ppat.1002882](https://doi.org/10.1371/journal.ppat.1002882) |
| **β-1,6-glucan** | DEFER (fold into β-1,3) | — | 0 | keep folded for v1 | — | — | Singh et al. 2025, *mSphere*, PMID [40920062](https://pubmed.ncbi.nlm.nih.gov/40920062/), [DOI 10.1128/msphere.00341-25](https://doi.org/10.1128/msphere.00341-25) |
| **Mannan** | REDUCE coefficient | `cpd11685_c0` | 0.808 | **0.20–0.40** | ➖ coefficient-only | (no new rxn) — existing: `rxn09511_r0` GDP-Man → mannan in ER | Yeast mannoproteins ~25% wall vs filamentous-fungus ~10%. Bowman & Free 2006 |
| **1,3-β-glucan** | KEEP | `cpd11791_c0` | 1.135 | **~1.0** | ➖ coefficient-only | (no new rxn) — existing: `rxn09400_c0` UDP-Glc → cpd00014 + cpd00067 + **cpd11791** | Major wall polysaccharide in both yeast and filamentous fungi |

## 2. Storage / osmolytes

| Compound | Action | Metabolite cpd | Current | Suggested | Pathway | MS rxn IDs + equations | Evidence |
|---|---|---|---|---|---|---|---|
| **D-Mannitol** | ADD to biomass + gap-fill 2 reactions | `cpd00314_c0` (mannitol) + `cpd27436_c0` (mannitol-1-P) | 0 | **0.10–0.50** (start ~0.20) | 🆕 2 rxns | ① **rxn00546** (EC 1.1.1.17) D-Mannitol-1-P:NAD⁺ 5-oxidoreductase: cpd00003 (NAD) + **cpd27436** (mannitol-1-P) ⇌ cpd00004 (NADH) + cpd00067 + cpd00072 (F6P)<br>② **rxn01560** (EC 3.1.3.22) D-Mannitol-1-P phosphohydrolase: cpd00001 + **cpd27436** ⇌ cpd00009 + **cpd00314** (mannitol) | Solomon, Tan & Oliver 2005, *MPMI* 18:110-115, PMID [15720079](https://pubmed.ncbi.nlm.nih.gov/15720079/), [DOI 10.1094/MPMI-18-0110](https://doi.org/10.1094/MPMI-18-0110). Compatible-solutes review: Dijksterhuis et al. 2006, *Biochem J*, PMID [16987106](https://pubmed.ncbi.nlm.nih.gov/16987106/), [DOI 10.1042/BJ20061229](https://doi.org/10.1042/BJ20061229) |
| **Trehalose** | KEEP | `cpd00794_c0` | 0.0234 | keep | ➖ | (no new rxn) — `rxn02004_c0` TPS/TPP ✅ in model | Universal fungal storage carb |
| **Glycogen** | KEEP | `cpd00155_c0` | 0.519 | keep ~0.3–0.5 | ➖ | (no new rxn) — `rxn13300_c0` UDP-Glc → glycogen ✅ in model | Yeast value reasonable for filamentous fungi |
| **Glycerol (appressorial)** | DO NOT ADD to vegetative biomass | `cpd00100_c0` | — | — | ➖ | (no new rxn) — biosynthesis ✅ in model | Howard et al. 1991, *PNAS* 88:11281-11284, PMID [1837147](https://pubmed.ncbi.nlm.nih.gov/1837147/), [DOI 10.1073/pnas.88.24.11281](https://doi.org/10.1073/pnas.88.24.11281). Qin et al. 2023, *Int J Mol Sci* 24:7411, PMID [37108573](https://pubmed.ncbi.nlm.nih.gov/37108573/) |

## 3. Membrane / lipids

| Compound | Action | Metabolite cpd | Current | Suggested | Pathway | MS rxn IDs + equations | Evidence |
|---|---|---|---|---|---|---|---|
| **Zymosterol** | REDUCE to ~0 | `cpd03221_c0` | 0.0015 | **≤1e-5** | ➖ coefficient-only | (no new rxn) — squalene→…→zymosterol pathway ✅ in model | Intermediate, not membrane sterol |
| **Ergosterol** | INCREASE (absorb zymosterol's share) | `cpd01170_c0` | 0.0001 | **0.0015–0.0020** | ➖ coefficient-only | (no new rxn) — full pathway to ergosterol ✅ in model | Sun et al. 2025, *J Agric Food Chem* 73:8290-8302, PMID [40298938](https://pubmed.ncbi.nlm.nih.gov/40298938/), [DOI 10.1021/acs.jafc.5c02145](https://doi.org/10.1021/acs.jafc.5c02145) |
| **Yeast MIPC / M(IP)₂C sphingolipids (7 species)** → **fungal β-glucosylceramide** | REMOVE 7 yeast species; ADD 1 fungal GlcCer | drop cpd15253, cpd15254, cpd15255, cpd15256, cpd15257, cpd15258, cpd15264; ADD `cpd00878_c0` (D-glucosylceramide; in ModelSEED, not in fsp237) | sum ≈ 0.012 | drop yeast species; **add cpd00878 at ~0.005–0.010** | 🆕 4 rxns (Δ4-desat, Δ8-desat, C9-MT, GCS) | ① **rxn44415** (EC 1.14.19.17) Δ4-sphinganine/dihydroceramide desaturase: cpd00007 (O₂) + 2 cpd27031 (cytochrome b5 red) + **cpd36512** (sphinganine acyl-CoA pool) ⇌ 2 cpd00001 + cpd00878 (Cer) + 2 cpd27029<br>② **rxn43251** (EC 1.14.19.29) Δ8-sphingolipid desaturase: cpd00007 + 2 cpd27031 + **cpd27876** ⇌ 2 cpd00001 + 2 cpd27029 + **cpd36798** (8-sphingenine analog)<br>③ **C9 methyltransferase** (EC 2.1.1.-): SAM (cpd00017) + 4,8-sphingadiene → SAH (cpd00019) + 9-Me-4,8-sphingadiene. *ModelSEED has no curated entry — add as ad-hoc rxn with new compound stub if you want the methyl branch explicit; otherwise lump into the GCS step.*<br>④ **rxn01088** (EC 2.4.1.80) UDP-glucose:N-acylsphingosine D-glucosyltransferase: cpd00026 (UDP-Glc) + cpd00167 (ceramide) ⇌ cpd00014 (UDP) + **cpd00878** (glucosylceramide) | Toledo et al. 1999, *Biochemistry* 38:7294-7306, PMID [10353841](https://pubmed.ncbi.nlm.nih.gov/10353841/), [DOI 10.1021/bi982898z](https://doi.org/10.1021/bi982898z). Follow-up: Toledo et al. 2000, *J Lipid Res*, PMID [10787440](https://pubmed.ncbi.nlm.nih.gov/10787440/) |
| **PC, PA, PS, PE, TAG** | KEEP | cpd11624, cpd15276, cpd29687, cpd29688, cpd11677 | each ~1e-4 | keep | ➖ | (no new rxn) — headgroup biosynthesis ✅ in model | Phospholipid headgroups conserved across fungi |

## 4. Cofactors / quinones

| Compound | Action | Metabolite cpd | Current | Suggested | Pathway | MS rxn IDs + equations | Evidence |
|---|---|---|---|---|---|---|---|
| **Ubiquinone-6 → Ubiquinone-9** | REPLACE | drop `cpd15290_m0` (Q6); ADD `cpd01351_m0` (Q9; in ModelSEED, not in fsp237) | 0.0002 (Q6) | **0.0002 (Q9)** | 🆕 3 rxns (chain extension + head condensation + final ring methylation) | ① **rxn16110** (EC 2.5.1.84) GPP:IPP transtransferase (adds 7 IPP): 7 cpd00113 (IPP) + cpd00283 (GPP) → 7 cpd00012 + 7 cpd00067 + **cpd02172** (all-trans-nonaprenyl/solanesyl-PP)<br>② **rxn05034** (EC 2.5.1.39) Solanesyl-PP:4-HB nonaprenyltransferase: cpd00136 (4-HB) + cpd02172 → cpd00012 + cpd00067 + **cpd02419** (3-nonaprenyl-4-HB)<br>③ **rxn05003** (EC 2.1.1.64) SAM:nonaprenyl ring methyltransferase (last Coq step): cpd00017 (SAM) + cpd02064 ⇌ cpd00019 (SAH) + **cpd01351** (Ubiquinone-9)<br>(Intermediate Coq2–Coq7 steps between ② and ③ can be lumped into one balancing reaction with O₂, NADH and SAM cofactors if you don't want to mirror the full Q6 cascade.) | Q-system survey: Kuraishi et al. 2000, *Antonie Van Leeuwenhoek* 77:179-186, PMID [10768477](https://pubmed.ncbi.nlm.nih.gov/10768477/), [DOI 10.1023/A:1002416431944](https://doi.org/10.1023/A:1002416431944). **Note:** no Colletotrichum-specific Q measurement found; Q9 is the most defensible interpolation from Sordariomycetes |
| **All other cofactors** | KEEP | NAD(P), CoA, FAD, PLP, SAM, THF, 5-Me-THF, siroheme, hemeA, riboflavin, TPP, biotin, pantothenate (14 cpds) | as-is | keep | ➖ | (no new rxn) — biosynthesis ✅ in model | Cofactors conserved across fungi |

## 5. Nucleotides — GC-content rebalance

| Component | Action | Metabolite cpds | Current (yeast ~38% GC) | Suggested (52.70% GC) | Pathway | MS rxn IDs + equations | Evidence |
|---|---|---|---|---|---|---|---|
| **dNTPs (dAMP, dTMP, dGMP, dCMP)** | REBALANCE | cpd00294, cpd00298, cpd00296, cpd00206 | dAMP 0.0036, dTMP 0.0036, dGMP 0.0024, dCMP 0.0024 (sum 0.012) | dAMP **0.00284**, dTMP **0.00284**, dGMP **0.00316**, dCMP **0.00316** | ➖ coefficient-only | (no new rxn) — ribonucleotide reductase + TS + kinases ✅ in model | Baroncelli et al. 2014, *Genome Announc* 2:e00540-14, PMID [24926053](https://pubmed.ncbi.nlm.nih.gov/24926053/), [DOI 10.1128/genomeA.00540-14](https://doi.org/10.1128/genomeA.00540-14). Sibling-species: Buiate et al. 2017, *BMC Genomics* 18:67, PMID [28073340](https://pubmed.ncbi.nlm.nih.gov/28073340/), [DOI 10.1186/s12864-016-3457-9](https://doi.org/10.1186/s12864-016-3457-9) |
| **NMPs (UMP, AMP, GMP, CMP)** | OPTIONAL rebalance | cpd00091, cpd00018, cpd00126, cpd00046 | UMP 0.0599, AMP 0.046, GMP 0.046, CMP 0.0447 | minor scaling optional | ➖ | (no new rxn) | RNA composition dominated by rRNA; genomic-GC effect small |

## 6. Amino acids — no change unless FSP237 proteomics is available

| Component | Action | Metabolite cpds | Current | Suggested | Pathway | MS rxn IDs + equations | Evidence |
|---|---|---|---|---|---|---|---|
| **20 standard amino acids** | KEEP | 20 cpds | yeast (Ala 0.4588 → Cys 0.0066) | keep | ➖ | (no new rxn) — biosynthesis ✅ in model | Avg protein AA composition varies <10% across fungi |

## 7. Energy / Growth-Associated Maintenance (GAM)

| Component | Action | Metabolite cpds | Current | Suggested | Pathway | MS rxn IDs + equations | Evidence |
|---|---|---|---|---|---|---|---|
| **ATP + H₂O → ADP + Pi + H⁺** | KEEP | cpd00001, cpd00002 → cpd00008, cpd00009, cpd00067 | 59.28 | keep ~59 | ➖ | (no new rxn) | Yeast GAM in same ballpark as filamentous fungi |

## 8. Pathogenicity-stage compounds (not in vegetative biomass on glucose)

| Compound | Action | Metabolite cpds | Current | Suggested | Pathway | MS rxn IDs + equations | Evidence |
|---|---|---|---|---|---|---|---|
| **DHN-melanin** | DEFER to appressorial-stage biomass | needs new polymer cpd + intermediates `cpd02495` (1,3,6,8-THN), `cpd00578` (scytalone), `cpd00864` (1,3,8-THN), `cpd01141` (vermelone), `cpd_DHN`, `cpd_melanin` | 0 | 0 for vegetative | 🆕 ~5 rxns (only if you build appressorial-stage biomass) | ① **rxn41997** (EC 2.3.1.233) THN synthase / PKS1: 5 cpd00067 + 5 cpd00070 (malonyl-CoA) ⇌ cpd00001 + 5 cpd00010 (CoA) + 5 cpd00011 (CO₂) + **cpd02495** (1,3,6,8-THN)<br>② **rxn02079** (EC 1.1.1.252) Scytalone:NADP⁺ Δ5-oxidoreductase: cpd00006 + **cpd00578** (scytalone) → cpd00005 + cpd00067 + cpd02495 (i.e. *runs in reverse* in melanin: 1,3,6,8-THN → scytalone via 4-HNR; same enzyme)<br>③ **rxn02080** (EC 4.2.1.94) Scytalone 7,8-hydro-lyase: cpd00578 → cpd00001 + **cpd00864** (1,3,8-THN)<br>④ **rxn02381** (EC 1.1.1.252) Vermelone:NADP⁺ Δ5-oxidoreductase: cpd00006 + **cpd01141** (vermelone) → cpd00005 + cpd00067 + cpd00864 (i.e. *runs in reverse*: 1,3,8-THN → vermelone)<br>⑤ Vermelone dehydratase → 1,8-DHN, then laccase polymerization → DHN-melanin. *No curated ModelSEED entries for the last two steps — add as ad-hoc reactions with new compound stubs.* | Qin et al. 2023, *Int J Mol Sci* 24:7411, PMID [37108573](https://pubmed.ncbi.nlm.nih.gov/37108573/), [DOI 10.3390/ijms24087411](https://doi.org/10.3390/ijms24087411). Melanin biosyn in *C. gloeosporioides*: Wang et al. 2021, *Fungal Biol* 125:813-824, PMID [34420695](https://pubmed.ncbi.nlm.nih.gov/34420695/), [DOI 10.1016/j.funbio.2021.04.004](https://doi.org/10.1016/j.funbio.2021.04.004) |

---

## Appendix — Compound ID legend for all proposed new metabolites

| cpd ID | Name | Status |
|---|---|---|
| cpd11683 | Chitin | ✅ already in fsp237 |
| cpd12148 | α-1,3-D-glucan (α-D-glucan polymer) | in ModelSEED, **needs adding** to fsp237 |
| cpd00314 | D-Mannitol | in ModelSEED, **needs adding** to fsp237 |
| cpd27436 | D-Mannitol-1-phosphate | in ModelSEED, **needs adding** to fsp237 |
| cpd00878 | D-Glucosylceramide (β-GlcCer) | in ModelSEED, **needs adding** to fsp237 |
| cpd00167 | N-acylsphingosine (ceramide substrate for GCS) | check if in fsp237 |
| cpd01351 | Ubiquinone-9 | in ModelSEED, **needs adding** to fsp237 |
| cpd02172 | All-trans-nonaprenyl-diphosphate (solanesyl-PP) | in ModelSEED, **needs adding** to fsp237 |
| cpd02419 | 3-Nonaprenyl-4-hydroxybenzoate | intermediate, **needs adding** if not lumping |
| cpd02064 | 2-Methoxy-6-nonaprenyl-1,4-benzoquinol (Q9 final intermediate) | intermediate, **needs adding** if not lumping |
| cpd02495 | 1,3,6,8-tetrahydroxynaphthalene (THN) | in ModelSEED |
| cpd00578 | Scytalone | in ModelSEED |
| cpd00864 | 1,3,8-trihydroxynaphthalene | in ModelSEED |
| cpd01141 | Vermelone | in ModelSEED |
| cpd_DHN, cpd_melanin | 1,8-DHN and DHN-melanin polymer | **needs ad-hoc creation** |

## Summary — gap-fill workload

| Compound | New ModelSEED reactions | New metabolites |
|---|---|---|
| Chitin | 0 (pathway complete) | 0 |
| α-1,3-glucan | 1 (rxn15561) | 1 (cpd12148) |
| Mannitol | 2 (rxn00546 + rxn01560) | 2 (cpd00314, cpd27436) |
| Glucosylceramide | 3 curated (rxn44415, rxn43251, rxn01088) + 1 ad-hoc (C9-MT) | ≥2 (cpd00878 + intermediates) |
| Ubiquinone-9 | 3 (rxn16110 + rxn05034 + rxn05003) — optionally lumped | 4 (cpd01351 + cpd02172 + cpd02419 + cpd02064) |
| Ergosterol/zymosterol | 0 | 0 |
| Nucleotide rebalance | 0 | 0 |
| Mannan reduction | 0 | 0 |
| DHN-melanin (deferred) | 4 curated (rxn41997, rxn02079, rxn02080, rxn02381) + 2 ad-hoc | ~6 |
| **TOTAL for v1** (chitin + mannitol + α-glucan + ergosterol swap + Q9 + GlcCer) | **9 curated + 1 ad-hoc** | **9** |

## What's still uncertain — flag before committing coefficients

1. **FSP237 vs TX430BB**: published reference is TX430BB. Treat the 52.70% GC and other genome stats as a proxy for FSP237.
2. **Ubiquinone form**: no *Colletotrichum*-specific measurement in PubMed; Q9 is the most defensible interpolation from Sordariomycetes/Hypocreales but should be checked against JGI MycoCosm metadata.
3. **α-1,3-glucan in *Colletotrichum* specifically**: AGS1 function firmly established for *Magnaporthe* (Fujikawa 2012); genomic evidence for AGS orthologs in *Colletotrichum* exists (Buiate 2017). No direct quantitative measurement in PubMed.
4. **Mannitol % dry weight in *Colletotrichum* mycelium**: pathway evidence + phenotypic role in *Stagonospora nodorum* (Solomon 2005); no direct *C. sublineola* number. Suggested 0.10–0.50 brackets published values for filamentous fungi (~2–10% dry weight).
5. **Fungal GlcCer fraction in *Colletotrichum* total lipid**: structural difference well-established (Toledo 1999); absolute quantitation not in literature pulled. Suggested coefficient anchored to current yeast-sphingolipid sum (0.012).
6. **C9-methyltransferase**: no curated ModelSEED entry for the fungal-specific sphingolipid C9-methyltransferase. If you want the methyl branch explicit, add as an ad-hoc reaction; otherwise the GCS step (`rxn01088`) lumps over it.
