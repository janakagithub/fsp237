# Peroxisomal β-oxidation (lipid / fatty-acid degradation)

## 1. Finding

Two independent signals point at peroxisomal β-oxidation in the FSP237 (*Colletotrichum
sublineola*) V10 GEM analysis:

- **Stage-differential (flux).** Active-rate is **0.1148 in pre-infection** media, **0.0 in
  biotrophic**, and **0.0 in necrotrophic** media (range 0.1148). The pathway carries flux
  **only** under pre-infection conditions.
- **Method-discordant.** It is the **2nd most discordant pathway** across integration
  methods (mean_method_range **0.4475**): pFBA, E-Flux, GIMME and iMAT disagree strongly on
  whether it is active.

Because there is a single static transcriptome (S1, pathogenic state) applied to every
medium, the stage-differential pattern is **flux/medium-driven, not expression-driven**
(see Caveats).

## 2. Mechanistic interpretation

The pre-infection-only flux is biologically the *expected* signature of an appressorium-
forming pathogen. Before penetration, the spore/germling has essentially no external carbon
and must run on **endogenous storage lipid** (triacylglycerol in lipid bodies). Mobilized
fatty acids are degraded by **peroxisomal β-oxidation** to acetyl-CoA, which is routed
through the **glyoxylate cycle** (isocitrate lyase / malate synthase) to make C4 units for
gluconeogenesis and, critically, the **glycerol** that generates the enormous appressorial
turgor used to breach the host cuticle. In hemibiotrophs this lipid economy is a
pre-penetration / early-infection program: once inside living host cells (biotrophy) and
later during necrotrophy, the fungus has access to host-derived sugars and its carbon source
shifts, so a medium encoding host nutrient availability correctly zeroes β-oxidation flux.
Thus "β-oxidation active only in pre-infection media" is consistent with the storage-lipid →
β-oxidation → glyoxylate cycle → turgor/gluconeogenesis axis that is repeatedly shown to be
required for appressorium function and virulence (Evidence below).

The method discordance is also expected rather than alarming. β-oxidation flux is highly
sensitive to (i) the **medium/carbon constraints** (whether any external carbon is provided,
and whether the glyoxylate-cycle bypass is open), and (ii) **which of several redundant
acyl-CoA-oxidase / multifunctional-protein genes** and their GPR (gene–protein–reaction)
rules the expression prior maps onto. With one flat expression vector, E-Flux (scales bounds
by expression), GIMME (penalizes low-expression reactions), iMAT (discretizes into on/off),
and unconstrained pFBA make different choices when the gene evidence is ambiguous and the
reaction is near a solution-space edge — producing a wide method range even though the
underlying biology (a switchable, medium-gated lipid catabolic route) is coherent.

## 3. Evidence

**Tier 1 — *Colletotrichum* spp. (anthracnose)**
- Kimura A. et al. 2001, *Plant Cell* — Peroxisomal metabolic function is required for
  appressorium-mediated plant infection by *Colletotrichum lagenarium* (= *C. orbiculare*).
  Establishes that peroxisome-based metabolism is essential for penetration. [Tier 1]
  DOI: 10.1105/tpc.010084
- Asakura M. et al. 2006, *Appl. Environ. Microbiol.* — Multiple contributions of
  peroxisomal metabolic function (fatty-acid β-oxidation among them) to fungal pathogenicity
  in *Colletotrichum lagenarium*. Directly ties β-oxidation/peroxisomal metabolism to
  virulence. [Tier 1] DOI: 10.1128/AEM.00988-06
- Fujihara N. et al. 2010, *Mol. Plant-Microbe Interact.* — Peroxisome-biogenesis factor
  **PEX13** required for appressorium-mediated infection by *C. orbiculare*; PEX mutants that
  cannot import peroxisomal β-oxidation enzymes are penetration-defective. [Tier 1]
  DOI: 10.1094/MPMI-23-4-0436
- Kubo Y. et al. 2015, *mBio* — *C. orbiculare* **FAM1** encodes a Woronin-body-associated
  **Pex22** peroxin required for appressorium-mediated infection. [Tier 1]
  DOI: 10.1128/mBio.01305-15
- Wang X. et al. 2026, *J. Fungi (Basel)* — **Pex8**, a fungal-specific peroxin, regulates
  peroxisome biogenesis and pathogenicity in the cucumber anthracnose fungus *C. orbiculare*.
  Recent confirmation that peroxisome function underpins *Colletotrichum* virulence. [Tier 1]
  DOI: 10.3390/jof12040248

**Tier 2 — related appressorial / hemibiotrophic pathogens (esp. *Magnaporthe oryzae*)**
- Wang Z.-Y. & Talbot N.J. 2007, *Mol. Plant-Microbe Interact.* — Functional analysis of
  lipid metabolism in *M. grisea* reveals a requirement for **peroxisomal fatty-acid
  β-oxidation** during appressorium-mediated plant infection. The most direct Tier-2 support
  for this pathway's role. [Tier 2] DOI: 10.1094/MPMI-20-5-0475
- Wang Z.-Y. & Talbot N.J. (Wang et al.) 2003, *Mol. Microbiol.* — The **glyoxylate cycle**
  is required for temporal regulation of virulence in *M. grisea*, linking β-oxidation-derived
  acetyl-CoA to infection. [Tier 2] DOI: 10.1046/j.1365-2958.2003.03412.x
- Thines E., Weber R.W.S. & Talbot N.J. 2000, *Plant Cell* — MAP-kinase/PKA-dependent
  mobilization of **triacylglycerol and glycogen** during appressorium turgor generation;
  storage-lipid mobilization is the upstream step feeding β-oxidation. [Tier 2]
  DOI: 10.1105/tpc.12.9.1703
- Bhadauria V. et al. 2012, *PLoS ONE* — Peroxisomal alanine:glyoxylate aminotransferase
  **AGT1** is indispensable for appressorium function of *M. oryzae* (triglyceride
  mobilization/utilization during infection). [Tier 2] DOI: 10.1371/journal.pone.0036266
- Chen X.-L. et al. 2017, *Mol. Plant Pathol.* — **Peroxisomal fission** is induced during
  appressorium formation and is required for full virulence of the rice blast fungus. [Tier 2]
  DOI: 10.1111/mpp.12395

**Tier 3 — reviews / general fungal pathway logic**
- Falter C. & Reumann S. 2022, *Mol. Plant Pathol.* — The essential role of fungal
  peroxisomes in plant infection (synthesizes β-oxidation/glyoxylate-cycle → appressorium
  turgor across taxa). [Tier 3] DOI: 10.1111/mpp.13180
- Kubo Y. 2013, *Subcell. Biochem.* — Function of peroxisomes in plant–pathogen
  interactions. [Tier 3] DOI: 10.1007/978-94-007-6889-5_18

## 4. Caveats

- **Single-transcriptome limitation.** Only S1 (pathogenic state, 3 technical reps) exists;
  it is applied unchanged to all 18 media. The stage-differential signal therefore reflects
  **medium/constraint differences**, not measured stage-specific expression. We cannot claim
  β-oxidation is *transcriptionally* induced pre-infection — only that the model routes flux
  there when the pre-infection medium (low external carbon) forces reliance on stored lipid.
  Biologically the literature supports pre-penetration lipid catabolism, but this dataset
  cannot prove it for FSP237.
- **Model/method artifacts.** The very method discordance (0.4475) means "activity" here is
  method-dependent. β-oxidation sits near a solution-space edge and depends on how each method
  treats ambiguous GPRs for the redundant acyl-CoA-oxidase / multifunctional-β-oxidation
  genes; flux may be an alternate-optimum choice (pFBA parsimony, iMAT on/off thresholding)
  rather than a robust prediction. Treat the 0.1148 magnitude as qualitative.
- **Alternatives.** A zero in biotrophic/necrotrophic media could equally reflect the medium
  supplying carbon that makes β-oxidation unnecessary, or the objective function preferring
  glycolytic carbon — not necessarily that the fungus shuts β-oxidation off in planta.
- **Taxonomic distance.** Strongest *Colletotrichum* β-oxidation evidence is from
  *C. lagenarium/orbiculare* (cucurbit anthracnose); *C. sublineola* is a graminicolous
  relative, so cross-species transfer is reasonable but not identical.

## 5. Confidence

**Medium–High.** The mechanistic story (pre-infection storage-lipid → peroxisomal
β-oxidation → glyoxylate cycle → appressorium turgor/virulence) is exceptionally well
supported by direct *Colletotrichum* genetics and by *Magnaporthe* β-oxidation mutants, so
the *direction* of the pre-infection signal is highly credible. Confidence is held below
"High" because the flux magnitude rests on a single static transcriptome and the pathway is
the 2nd-most method-discordant, making the quantitative flux value itself unreliable.

```json
{
  "pathway": "Peroxisomal β-oxidation (lipid/fatty-acid degradation)",
  "one_liner": "Pre-infection-only flux matches the storage-lipid → peroxisomal β-oxidation → glyoxylate-cycle → appressorium-turgor program required for penetration in Colletotrichum and Magnaporthe; strong mechanistic support but flux magnitude is uncertain (single static transcriptome, 2nd-most method-discordant).",
  "confidence": "Medium-High",
  "top_refs": [
    {"cite": "Wang Z.-Y. & Talbot N.J. 2007, Mol. Plant-Microbe Interact. — peroxisomal fatty-acid β-oxidation required for appressorium-mediated infection in Magnaporthe grisea", "url": "https://doi.org/10.1094/MPMI-20-5-0475"},
    {"cite": "Kimura A. et al. 2001, Plant Cell — peroxisomal metabolic function required for appressorium-mediated infection by Colletotrichum lagenarium", "url": "https://doi.org/10.1105/tpc.010084"},
    {"cite": "Asakura M. et al. 2006, Appl. Environ. Microbiol. — multiple contributions of peroxisomal metabolic function to pathogenicity in Colletotrichum lagenarium", "url": "https://doi.org/10.1128/AEM.00988-06"},
    {"cite": "Wang Z.-Y. et al. 2003, Mol. Microbiol. — glyoxylate cycle required for temporal regulation of virulence in Magnaporthe grisea", "url": "https://doi.org/10.1046/j.1365-2958.2003.03412.x"},
    {"cite": "Falter C. & Reumann S. 2022, Mol. Plant Pathol. — the essential role of fungal peroxisomes in plant infection (review)", "url": "https://doi.org/10.1111/mpp.13180"}
  ]
}
```
