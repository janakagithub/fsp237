# Central Carbon Metabolism — Glycolysis/Gluconeogenesis, TCA Cycle, Pentose Phosphate Pathway

*Pathway cluster interpretation for the FSP237 (*Colletotrichum sublineola*, sorghum anthracnose) V10 GEM infection-stage flux analysis. Cluster theme: the metabolic backbone — uniformly **active** across stages, but the region where transcriptomics-constrained integration methods most **disagree**.*

## 1. Finding

Across the 18-condition media panel, central carbon metabolism is uniformly active at every
infection stage, yet it is also where the four flux methods (pFBA, E-Flux, GIMME, iMAT)
diverge most. Reported as *active-rate (pre-infection / biotrophic / necrotrophic ; method-discordance = mean_method_range)*:

| Subsystem | Pre-infection | Biotrophic | Necrotrophic | Discordance |
|---|---|---|---|---|
| **Pentose phosphate pathway (PPP)** | 0.64 | 0.61 | 0.65 | **0.3167** |
| **TCA cycle** | 0.4294 | 0.4471 | 0.4118 | **0.4215** (3rd most discordant subsystem) |
| **Glycolysis / gluconeogenesis** | 0.4276 | 0.4414 | 0.4448 | **0.3372** |

Two signals stand out. (i) **Stage stability**: active-rates barely move across pre-infection →
biotrophic → necrotrophic (e.g., PPP 0.61–0.65; glycolysis 0.4276–0.4448) — central carbon
stays "on" regardless of infection stage. (ii) **Method discordance**: these hubs top the
method-disagreement ranking, with the TCA cycle the 3rd most discordant subsystem in the whole
model. The PPP carries a notably high active baseline (~0.6) that the methods still disagree about.

## 2. Mechanistic interpretation

**Biological reading — why central carbon is "always on."** Glycolysis/gluconeogenesis, the TCA
cycle, and the PPP form the shared backbone that supplies energy (ATP), redox cofactors
(NADH, NADPH), and the 12 precursor metabolites for essentially all biomass. Any viable,
growing cell must run them, so uniform high activity across stages is exactly what is expected
for a pathogen that is metabolically active in every phase modeled. The near-flat stage profile
here is consistent with the analysis design: a single static pathogenic-state transcriptome is
applied to every medium, so backbone activity tracks the demand for growth/maintenance rather
than any stage-specific transcriptional switch (see Caveats).

The infection-relevant nuances are carried in *how* the backbone is used:

- **PPP → NADPH for biosynthesis and oxidative-stress defense.** The oxidative PPP is the primary
  cytosolic source of NADPH, needed both for anabolism (lipids, amino acids, nucleotides during
  invasive growth) and to regenerate glutathione/thioredoxin pools that detoxify the host
  oxidative burst (ROS) encountered during penetration and biotrophic colonization. A high,
  stable PPP active-rate (~0.6, the highest of the three subsystems) fits a pathogen that must
  keep NADPH-dependent redox buffering and biosynthesis running throughout infection. In the rice
  blast fungus, glucose-6-phosphate flux partitioning into the PPP is a genetically enforced
  checkpoint for infection: Tps1 senses G6P and, via the PPP/NADPH balance, licenses pathogenic
  development, and a downstream transketolase (a non-oxidative PPP enzyme) checkpoint governs
  biotrophic growth inside living host cells.

- **Gluconeogenesis for growth on non-sugar host carbon.** Host tissue is not a free sugar buffet;
  early infection structures and biotrophic hyphae often subsist on lipids, and on organic
  acids/amino acids, requiring gluconeogenesis (and, from lipid-derived acetyl-CoA, the glyoxylate
  bypass) to build hexose-phosphate for cell-wall and PPP precursors. The glycolysis/gluconeogenesis
  subsystem being active — and slightly rising toward the necrotrophic stage in this analysis — is
  consistent with flexible carbon sourcing as the medium (and thus available carbon) changes. This
  reversible glycolytic/gluconeogenic axis is precisely the kind of bidirectional route that
  inflates model flexibility.

- **TCA cycle as the flexible energy/redox and biosynthesis hub.** The TCA cycle supplies NADH for
  respiration and carbon skeletons (α-ketoglutarate, oxaloacetate) for amino-acid biosynthesis, and
  can run in truncated/branched modes depending on carbon source. Its position as the most
  discordant of the three (and 3rd overall) reflects its high connectivity and the availability of
  anaplerotic/cataplerotic and glyoxylate alternatives.

**Methodological reading — why central carbon is a stress-test for integration methods.** These
subsystems are (a) the most highly connected nodes in the network, (b) rich in *reversible* and
*alternative* routes (glycolysis vs. gluconeogenesis; oxidative vs. non-oxidative PPP; TCA vs.
glyoxylate; multiple isozymes), and (c) largely degenerate in flux space — many alternate optima
carry the same objective value. Transcriptomic-integration methods differ in *how* they turn one
gene-expression vector into constraints: E-Flux scales bounds continuously by expression; GIMME
penalizes flux through below-threshold reactions while protecting the objective; iMAT maximizes
the number of high-expression reactions carrying flux (and minimizes low-expression ones) as a
discrete optimization. In a highly connected, alternative-route-rich region, small changes in the
expression threshold flip many reactions on/off and reroute flux among equivalent paths, so the
methods land on different-but-equally-optimal solutions. This is exactly the behavior benchmark
studies report: method choice strongly affects which reactions are called active, methods often
disagree with each other, and no single method is universally best. The **iMAT active-count
inflation caveat** is directly relevant here: because iMAT's objective *rewards activating*
high-expression reactions, it tends to report more reactions carrying (often small) flux, which
mechanically raises active-rates and widens the mean_method_range specifically in dense hubs like
central carbon. In short, the discordance is largely a property of the network topology and the
methods' differing treatment of degeneracy — not evidence of a real stage-specific biological
switch.

## 3. Evidence

- **[Tier 1]** *Colletotrichum* lifestyle transitions are accompanied by large-scale metabolic/
  transcriptional reprogramming, establishing that central metabolism is remodeled across biotrophy
  → necrotrophy in this genus — O'Connell RJ *et al.* (2012), *Nature Genetics*, "Lifestyle
  transitions in plant pathogenic *Colletotrichum* fungi deciphered by genome and transcriptome
  analyses." DOI: 10.1038/ng.2372 (PMID 22885923).
- **[Tier 1]** In the closely related maize anthracnose pathogen *Colletotrichum graminicola*
  (sister species to sorghum-infecting *C. sublineola*), establishment of biotrophy is marked by
  distinct early transcriptional/metabolic events, supporting stage-linked but continuously active
  central metabolism — Torres MF *et al.* (2016), *BMC Genomics*, "A *Colletotrichum graminicola*
  mutant deficient in the establishment of biotrophy reveals early transcriptional events in the
  maize anthracnose disease interaction." DOI: 10.1186/s12864-016-2546-0 (PMID 26956617).
- **[Tier 2]** PPP/NADPH partitioning is a genetically enforced checkpoint for fungal infection:
  Tps1 senses glucose-6-phosphate and controls the PPP, NADPH balance, and virulence in
  *Magnaporthe oryzae* — Wilson RA *et al.* (2007), *EMBO Journal*, "Tps1 regulates the pentose
  phosphate pathway, nitrogen metabolism and fungal virulence." DOI: 10.1038/sj.emboj.7601795
  (PMID 17641690).
- **[Tier 2]** An NADPH-dependent genetic switch (downstream of G6P/PPP flux) gates plant infection,
  tying redox cofactor supply from central carbon to pathogenicity — Wilson RA *et al.* (2010),
  *PNAS*, "An NADPH-dependent genetic switch regulates plant infection by the rice blast fungus."
  DOI: 10.1073/pnas.1006839107 (PMID 21115813).
- **[Tier 2]** A non-oxidative PPP enzyme (transketolase) checkpoint specifically governs biotrophic
  growth inside living rice cells, linking PPP flux to the biotrophic phase — Fernandez J,
  Marroquin-Guzman M, Wilson RA (2014), *PLoS Pathogens*, "Evidence for a transketolase-mediated
  metabolic checkpoint governing biotrophic growth in rice cells by the blast fungus *Magnaporthe
  oryzae*." DOI: 10.1371/journal.ppat.1004354 (PMID 25188286).
- **[Tier 2]** Review of central-carbon/nutrient strategy during hemibiotrophic invasion (glucose
  sensing, gluconeogenesis from lipid/organic-acid carbon, PPP/NADPH for redox) — Fernandez J,
  Wilson RA (2018), *Trends in Microbiology*, "Rise of a Cereal Killer: The Biology of *Magnaporthe
  oryzae* Biotrophic Growth." DOI: 10.1016/j.tim.2017.12.007 (PMID 29395728).
- **[Tier 2/3]** Biotrophic vs. necrotrophic lifestyles show niche-specific metabolic adaptation
  and differential nutrient/carbon use (filamentous plant pathogen, oomycete comparator) — Ah-Fong
  AMV *et al.* (2019), *PLoS Pathogens*, "Niche-specific metabolic adaptation in biotrophic and
  necrotrophic oomycetes is manifested in differential use of nutrients." DOI:
  10.1371/journal.ppat.1007729 (PMID 31002734).
- **[Methods]** Systematic benchmark showing transcriptomic-integration methods (incl. E-Flux,
  GIMME, iMAT) frequently disagree and none is universally superior — supports the discordance
  interpretation — Machado D, Herrgård M (2014), *PLoS Computational Biology*, "Systematic
  evaluation of methods for integration of transcriptomic data into constraint-based models of
  metabolism." DOI: 10.1371/journal.pcbi.1003580 (PMID 24762745).
- **[Methods]** Systematic evaluation showing that context-specific model outputs depend heavily on
  the extraction method and thresholds, especially for highly connected reactions — Opdam S *et al.*
  (2017), *Cell Systems*, "A Systematic Evaluation of Methods for Tailoring Genome-Scale Metabolic
  Models." DOI: 10.1016/j.cels.2017.01.010 (PMID 28215528).
- **[Methods]** Review of metabolic-modelling approaches for plant–microbe interactions, framing
  how expression data are integrated and their limitations — Feierabend M *et al.* (2025), *FEMS
  Microbiology Reviews*, "In silico encounters: harnessing metabolic modelling to understand
  plant–microbe interactions." DOI: 10.1093/femsre/fuaf030 (PMID 40705360).

## 4. Caveats

- **Single-transcriptome limitation (primary caveat).** Only one transcriptome (S1, pathogenic
  state, 3 technical replicates) is applied as a *static* expression prior to every medium.
  Expression does not vary by stage; only the medium (and hence flux) does. Therefore the
  near-flat stage profiles of central carbon are *flux/medium-driven*, not evidence of a
  transcriptional stage switch, and no differential-expression p-values exist for these numbers.
  The stage columns should be read as "how the backbone responds when the same expression prior
  meets different media," not as biological stage transitions.
- **Discordance is partly an artifact of network topology + method design.** The high mean_method_range
  in central carbon reflects reversibility, alternative optima, and hub connectivity as much as any
  biology. Different methods select different equally-optimal solutions; this is expected and is not
  a biological finding on its own.
- **iMAT active-count inflation.** iMAT's objective rewards activating high-expression reactions,
  which mechanically inflates the count of reactions carrying (possibly trivial) flux and widens the
  method range specifically in dense subsystems — so the PPP's high active-rate (~0.6) and the TCA
  discordance may be partly method-inflated. Active-rate should be interpreted as "reaction called
  carrying nonzero flux," not flux magnitude.
- **Directionality ambiguity.** Glycolysis and gluconeogenesis share reactions/annotations; an
  "active" call in this subsystem does not by itself indicate net glycolytic vs. gluconeogenic
  direction. Distinguishing them requires inspecting flux sign under each medium, not the active-rate.
- **Alternatives not run.** Flux variability analysis (FVA), flux sampling, or a consensus/ensemble
  across methods would quantify how much of the discordance is alternate-optima degeneracy vs. real
  method divergence; a stage-resolved transcriptome would be needed to test genuine stage effects.

## 5. Confidence

**Medium–High.** The *methodological* interpretation (central carbon as a discordance hotspot due
to connectivity, reversibility, and method design, incl. iMAT inflation) is High confidence —
strongly supported by benchmark literature and consistent with the reported numbers. The
*biological* interpretation (PPP/NADPH for redox defense and biosynthesis; gluconeogenesis on
host-derived non-sugar carbon; TCA flexibility) is Medium confidence — well supported in related
hemibiotrophs (*Magnaporthe*) and in *Colletotrichum* genus-level reprogramming, but not directly
demonstrated in *C. sublineola*, and dampened by the single-static-transcriptome design, which
prevents any stage-differential claim.

```json
{
  "pathway": "Central carbon metabolism (glycolysis/gluconeogenesis, TCA cycle, pentose phosphate pathway)",
  "one_liner": "Central carbon is uniformly active across all infection stages (PPP notably high, ~0.6) but is the top method-discordance region (TCA 3rd most discordant overall) because these highly connected, reversible, alternative-route hubs make transcriptomic-integration methods pick different equally-optimal solutions; biologically it supplies NADPH for redox/oxidative-stress defense and biosynthesis and enables gluconeogenesis on host non-sugar carbon, while the discordance and PPP active-rate are partly artifacts of network topology and iMAT active-count inflation under a single static transcriptome.",
  "confidence": "Medium-High",
  "top_refs": [
    {"cite": "Machado D, Herrgard M (2014) PLoS Comput Biol — Systematic evaluation of methods for integration of transcriptomic data into constraint-based models of metabolism", "url": "https://doi.org/10.1371/journal.pcbi.1003580"},
    {"cite": "Wilson RA et al. (2007) EMBO J — Tps1 regulates the pentose phosphate pathway, nitrogen metabolism and fungal virulence", "url": "https://doi.org/10.1038/sj.emboj.7601795"},
    {"cite": "Fernandez J, Marroquin-Guzman M, Wilson RA (2014) PLoS Pathog — Transketolase-mediated metabolic checkpoint governing biotrophic growth in rice cells", "url": "https://doi.org/10.1371/journal.ppat.1004354"},
    {"cite": "O'Connell RJ et al. (2012) Nat Genet — Lifestyle transitions in plant pathogenic Colletotrichum fungi deciphered by genome and transcriptome analyses", "url": "https://doi.org/10.1038/ng.2372"},
    {"cite": "Opdam S et al. (2017) Cell Syst — A Systematic Evaluation of Methods for Tailoring Genome-Scale Metabolic Models", "url": "https://doi.org/10.1016/j.cels.2017.01.010"}
  ]
}
```
