# Plant cell-wall degradation (necrotroph-specific)

*Pathway cluster: Pectin / D-galacturonate catabolism (Ashwell/uronate route) + pentose-sugar catabolism (hemicellulose-derived L-arabinose / D-xylose). FSP237 = *Colletotrichum sublineola*, sorghum anthracnose (hemibiotroph).*

## 1. Finding

Two related plant-cell-wall (PCW) catabolic buckets are predicted active **only** in necrotrophic-stage media, and are the sharpest stage-differential signals in the panel:

| Pathway | pre-infection | biotrophic | necrotrophic | range |
|---|---|---|---|---|
| Pectin / D-galacturonate catabolism ("Ashwell" pathway) | 0.0 | 0.0 | 0.1667 | 0.1667 |
| Pentose-sugar catabolism (L-arabinose / D-xylose) | 0.0 | 0.0 | 0.2 | 0.2 |

Pentose catabolism is the single most stage-differential pathway in the whole analysis (range 0.2). Both go from fully off (pre-infection, biotrophic) to on in necrotrophic media.

## 2. Mechanistic interpretation

Because the analysis applies ONE static expression prior (transcriptome S1) to every medium, these on/off differences are **substrate-availability (medium/flux) driven, not expression-driven**. The necrotrophic media evidently supply pectin- and hemicellulose-derived monomers — D-galacturonate (the dominant monomer of homogalacturonan pectin) and the pentoses L-arabinose and D-xylose (from arabinoxylan/arabinan hemicellulose) — while the pre-infection and biotrophic media do not. Given the catabolic genes are "expressed" in the static prior, the model simply routes flux through them whenever the corresponding substrate is present in the medium.

Biologically this is exactly the expected metabolic signature of the anthracnose **necrotrophic switch**. During early biotrophy *Colletotrichum* keeps host cells alive, deploys bulbous intracellular hyphae, and conspicuously *withholds* lytic plant-cell-wall-degrading enzymes (PCWDEs); at the transition to necrotrophy it releases a coordinated wave of pectinases (polygalacturonases, pectin/pectate lyases, pectin methylesterases) plus hemicellulases and cellulases that macerate host tissue, producing the water-soaked, collapsing anthracnose lesion. Cell-wall breakdown liberates D-galacturonate and pentoses, which the fungus then catabolizes as carbon/energy sources — the intracellular metabolic arm downstream of maceration. The flux model therefore recapitulates the correct *coupling*: necrotroph-stage substrate liberation → uronate + pentose catabolism switched on.

A biochemical caveat on the "Ashwell" label (see Caveats): filamentous fungi almost universally catabolize D-galacturonate through the **reductive** pathway (D-galacturonate reductase → L-galactonate → L-galactonate dehydratase → 2-keto-3-deoxy-L-galactonate → aldolase → pyruvate + L-glyceraldehyde; genes gaaA–gaaD), NOT the bacterial **oxidative** isomerase route of Ashwell / De Ley–Doudoroff (uronate isomerase → tagaturonate/altronate → KDG → KDG kinase → KDPG aldolase). Pentoses feed the fungal pentose catabolic pathway (L-arabinose/D-xylose reductases → xylitol/arabitol → D-xylulose → D-xylulose-5-P → pentose phosphate pathway). The direction of the switch is the robust result; the exact reaction naming is a model-annotation detail.

## 3. Evidence

- **[Tier 1]** *Colletotrichum* genomes are enriched for and temporally program PCWDE/CAZyme deployment across the biotrophy→necrotrophy switch, with pectinases and other lytic enzymes peaking at the necrotrophic transition — the landmark genome/transcriptome study of *C. graminicola* (maize) and *C. higginsianum*. O'Connell et al. 2012, *Nat Genet*. DOI 10.1038/ng.2372 (PMID 22885923).
- **[Tier 1]** Directly the study organism's sibling pair: comparative genomics of the host-specific siblings *C. graminicola* (maize) and *C. sublineola* (sorghum) catalogs their pathogenicity-gene and carbohydrate-active-enzyme / cell-wall-degrading repertoires. Buiate et al. 2017, *BMC Genomics*. DOI 10.1186/s12864-016-3457-9 (PMID 28073340).
- **[Tier 1]** Directly *C. sublineola*–sorghum: dual RNA-seq of the interaction captures in planta fungal gene expression, including secreted/degradative functions across infection. Vela et al. 2024, *Front Fungal Biol*. DOI 10.3389/ffunb.2024.1437344 (PMID 39220294).
- **[Tier 1]** Review of the sorghum–*C. sublineola* interaction and the hemibiotrophic lifestyle framing host-tissue maceration in the necrotrophic phase. Abreha et al. 2021, *Front Plant Sci*. DOI 10.3389/fpls.2021.641969 (PMID 33959139).
- **[Tier 1]** In another *Colletotrichum* (C. falcatum, sugarcane red rot), CAZy analysis links pectinolytic and cellulolytic enzyme deployment to pathogenesis — genus-level support that pectinase induction accompanies tissue destruction. Prasanth et al. 2022, *3 Biotech*. DOI 10.1007/s13205-022-03113-6 (PMID 35127303).
- **[Tier 2/3]** The fungal D-galacturonate catabolic route is the **reductive** pathway (distinct from the bacterial oxidative Ashwell/isomerase route); definitive comparative review. Richard & Hilditch 2009, *Appl Microbiol Biotechnol*. DOI 10.1007/s00253-009-1870-6 (PMID 19159926).
- **[Tier 3]** The four-gene fungal reductive D-galacturonate pathway (gaaA–gaaD) characterized/reconstituted from *Aspergillus niger*. Biz et al. 2016, *Microb Cell Fact*, DOI 10.1186/s12934-016-0544-1 (PMID 27538689); Kuivanen et al. 2015, *Microb Cell Fact*, DOI 10.1186/s12934-014-0184-2 (PMID 25566698).
- **[Tier 3]** Fungal flexibility and pathway wiring for catabolizing plant-biomass-derived monomers, including galacturonate and pentoses. Chroumpi et al. 2021, *Front Bioeng Biotechnol*. DOI 10.3389/fbioe.2021.644216 (PMID 33763411).
- **[Tier 3]** Regulation of the fungal pentose catabolic pathway (L-arabinose / D-xylose utilization) in *A. niger*. Peng et al. 2025, *Curr Res Microb Sci*. DOI 10.1016/j.crmicr.2025.100482 (PMID 41127560).
- **[Tier 3]** Efficient reductive D-galacturonate metabolism characterized genome-wide/enzymatically in a non-ascomycete fungus (Rhodosporidium toruloides), corroborating pathway generality. Protzko et al. 2019, *mSystems*. DOI 10.1128/msystems.00389-19 (PMID 31848309).

## 4. Caveats

- **Single static transcriptome (S1).** Expression does not vary by stage; the on/off pattern is driven entirely by which monomers each medium supplies. This is a *consistency* result (the model does the biologically sensible thing when given necrotroph substrates), not evidence that *C. sublineola* up-regulates these genes at the necrotrophic switch — although Tier 1 literature (O'Connell 2012) independently supports that it does.
- **"Ashwell pathway" is likely a mis-/loose annotation.** The Ashwell (uronate isomerase / De Ley–Doudoroff, oxidative) route is the *bacterial* pathway; fungi use the reductive galacturonate pathway. If the V10 GEM literally encodes the bacterial isomerase reactions, the *flux* conclusion (galacturonate catabolism on in necrotrophy) still holds, but the reaction stoichiometry/gene mapping should be checked against fungal gaaA–gaaD orthologs in FSP237.
- **Small absolute active-rates (0.1667, 0.2)** reflect a minority of media/methods turning the pathway on; magnitudes are fractions of conditions, not flux quantities, so should not be over-interpreted.
- **Alternative explanation:** the signal could partly be a medium-composition artifact (necrotroph media were formulated with pectin/hemicellulose hydrolysates), rather than a de novo metabolic prediction — i.e., the model is partly reporting the media design back. The biological plausibility (matching the known maceration phase) is what makes it interpretable.
- Could not additionally verify a *C. graminicola*-specific PCWDE-wave transcriptome paper (repeated web queries hit an intermittent content filter); O'Connell 2012 and the two *C. sublineola* papers already cover the Tier 1 claim, so no fabricated reference was substituted.

## 5. Confidence

**Medium–High.** The direction of the result (pectin/D-galacturonate + pentose catabolism switched on specifically in necrotrophic media) is strongly consistent with well-established, *Colletotrichum*-specific biology of the necrotrophic maceration phase (Tier 1: O'Connell 2012, Buiate 2017, Vela 2024). Confidence is held below "High" only because the signal is substrate-driven under a single static transcriptome and because the "Ashwell" annotation may not match the fungal reductive pathway biochemistry.

```json
{
  "pathway": "Plant cell-wall degradation (pectin/D-galacturonate + pentose catabolism), necrotroph-specific",
  "one_liner": "Pectin/D-galacturonate and hemicellulose-pentose catabolism switch on only in necrotrophic media (range 0.1667 and 0.2), a substrate-driven signature matching the Colletotrichum necrotrophic tissue-maceration phase; note fungi use the reductive galacturonate pathway, not the bacterial oxidative Ashwell route.",
  "confidence": "Medium-High",
  "top_refs": [
    {"cite": "O'Connell et al. 2012, Nat Genet — Lifestyle transitions in plant pathogenic Colletotrichum (PCWDE/CAZyme waves at necrotrophic switch)", "url": "https://doi.org/10.1038/ng.2372"},
    {"cite": "Buiate et al. 2017, BMC Genomics — Comparative genomics of C. graminicola and C. sublineola pathogenicity genes/CAZymes", "url": "https://doi.org/10.1186/s12864-016-3457-9"},
    {"cite": "Vela et al. 2024, Front Fungal Biol — Dual RNA-seq of sorghum-C. sublineola interaction", "url": "https://doi.org/10.3389/ffunb.2024.1437344"},
    {"cite": "Richard & Hilditch 2009, Appl Microbiol Biotechnol — D-galacturonic acid catabolism in microorganisms (reductive fungal vs oxidative Ashwell route)", "url": "https://doi.org/10.1007/s00253-009-1870-6"}
  ]
}
```
