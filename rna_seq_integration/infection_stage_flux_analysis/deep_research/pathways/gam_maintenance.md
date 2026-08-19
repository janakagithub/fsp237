# Growth-Associated Maintenance (GAM) / Non-Growth Maintenance Energy

## 1. Finding

The **GAM / maintenance** category is the single **most method-discordant pathway** in the
FSP237 panel: **mean_method_range = 1.0** — the maximum possible value on the active-rate
scale [0, 1]. Inspection of `pathway_analysis/pathway_method_matrix.tsv` shows *why* this
number is maximal, and the pattern is diagnostic:

- The category is a **single reaction** (`n_rxn = 1`) — the ATP-maintenance / GAM-coupled
  ATP-hydrolysis reaction (of the general form `ATP + H2O -> ADP + Pi + H+`).
- **iMAT** carries a small flux (`total_flux ≈ 0.001`, `active_rate = 1.0`) in **all 18/18
  conditions**, across every medium, O2 level, and infection stage.
- **pFBA, E-Flux, and GIMME** carry **zero flux** (`active_rate = 0.0`) in **all 18/18
  conditions**.

So the "discordance" is a clean **binary split between methods** (iMAT ON everywhere vs. the
other three OFF everywhere), invariant to medium and stage. `max(active_rate) −
min(active_rate) = 1.0 − 0.0 = 1.0`. This is the maximal-disagreement signature, and it is
produced by the *methods*, not by the *biology* or even the *medium*.

## 2. Mechanistic interpretation

**This is a modeling artifact, and its cause is identifiable from the numbers.** Two model
features combine to force maximal disagreement:

1. **No gene–protein–reaction (GPR) rule.** ATP-maintenance and biomass-coupled GAM terms
   are *modeling constructs*, not enzyme-catalyzed reactions. They are added by convention to
   represent the aggregate energetic cost of staying alive (protein turnover, ion homeostasis,
   macromolecular repair) and of polymerizing biomass, and are lumped into a stoichiometric
   ATP-hydrolysis pseudo-reaction with **no associated gene** (Thiele & Palsson 2010; Feist &
   Palsson 2010; Varma & Palsson 1994). Every transcriptomics-integration method
   (E-Flux, GIMME, iMAT) maps expression onto reactions **through the GPR**. A reaction with
   no GPR receives **no expression signal**, so these methods have *nothing to constrain it
   with* and fall back on their default handling — which differs by method.

2. **The split tracks flux-parsimony, not expression.** With the maintenance reaction's
   lower bound at (or near) zero and no expression prior:
   - **pFBA** explicitly minimizes total flux after fixing the objective (Lewis et al. 2010),
     so a non-forced ATP-hydrolysis reaction is driven to **0**.
   - **E-Flux** bounds fluxes by expression; a no-GPR reaction gets a default bound and,
     with the model minimizing/limiting flux, resolves to **0** (Colijn et al. 2009).
   - **GIMME** minimizes flux through below-threshold reactions subject to a required
     metabolic objective; an unannotated maintenance reaction is not required and is driven to
     **0** (Becker & Palsson 2008).
   - **iMAT** solves a MILP that **maximizes the number of reactions whose ON/OFF state
     agrees with discretized expression, with no flux-minimization step** (Shlomi et al. 2008;
     Zur et al. 2010). An unclassified no-GPR reaction is left free, and the MILP returns an
     **arbitrary small nonzero flux (≈0.001)** — a numerically incidental value, not a
     biological rate.

The tell-tale that this is non-biological: the flux is **identical (0.001) across all 18
media and both O2 states** in iMAT and **identically zero** in the other three. A genuine
energetic signal would vary with carbon source, O2 availability, and infection stage. Here it
does not vary at all — the outcome is fully determined by each method's objective structure.
Because the category contains only one reaction, the pathway-level range collapses to the
reaction's binary ON/OFF disagreement, guaranteeing the extremal value of 1.0.

**Genuine biology, as context (not the headline).** Infection *is* energetically expensive,
and a properly parameterized maintenance/GAM term is biologically meaningful:
- Appressorium-mediated penetration in related fungi generates enormous turgor (up to ~8 MPa
  in *Magnaporthe oryzae*), built from glycerol accumulation — an osmotic/ATP-demanding
  process (Howard et al. 1991; de Jong et al. 1997).
- The hemibiotrophy-to-necrotrophy switch in *Colletotrichum* involves massive stage-specific
  secretion of effectors and CAZymes (O'Connell et al. 2012), and secretion, membrane
  transport, and protein turnover all draw on maintenance ATP.

But none of this can be read out of the GAM discordance here: with a static single
transcriptome and a no-GPR maintenance reaction, the model has **no channel** to express
stage-specific energetics through this category. The discordance is expected by construction.

## 3. Evidence

- **[Tier 3]** ATP-maintenance (NGAM) and growth-associated maintenance (GAM) are non-gene-
  associated pseudo-reactions added by reconstruction convention; the protocol prescribes
  estimating them from chemostat/growth data and encoding them as ATP-hydrolysis terms —
  Thiele I, Palsson BØ (2010), *A protocol for generating a high-quality genome-scale
  metabolic reconstruction*, Nature Protocols 5(1):93–121. DOI: 10.1038/nprot.2009.203
- **[Tier 3]** The biomass objective and its coupled GAM/ATP-maintenance terms are modeling
  abstractions, not enzymatic reactions — Feist AM, Palsson BO (2010), *The biomass objective
  function*, Curr Opin Microbiol 13(3):344–349. DOI: 10.1016/j.mib.2010.03.003
- **[Tier 3]** Original empirical GAM/NGAM estimation from growth data (the convention's
  basis) — Varma A, Palsson BO (1994), *Stoichiometric flux balance models quantitatively
  predict growth and metabolic by-product secretion in wild-type Escherichia coli W3110*,
  Appl Environ Microbiol 60(10):3724–3731. DOI: 10.1128/aem.60.10.3724-3731.1994
- **[Tier 3]** FBA/GPR framework — expression is mapped to flux *through* GPR rules, so
  no-GPR reactions receive no expression constraint — Orth JD, Thiele I, Palsson BØ (2010),
  *What is flux balance analysis?*, Nat Biotechnol 28(3):245–248. DOI: 10.1038/nbt.1614
- **[Tier 3]** Expression-integration methods disagree substantially and are sensitive to
  method choice/parameterization, especially for reactions lacking direct expression mapping —
  Machado D, Herrgård M (2014), *Systematic Evaluation of Methods for Integration of
  Transcriptomic Data into Constraint-Based Models of Metabolism*, PLoS Comput Biol
  10(4):e1003580. DOI: 10.1371/journal.pcbi.1003580
- **[Tier 3]** Independent systematic evaluation confirming method-dependent, often
  low-concordance outputs from context-specific model extraction — Opdam S et al. (2017),
  *A Systematic Evaluation of Methods for Tailoring Genome-Scale Metabolic Models*, Cell Syst
  4(3):318–329. DOI: 10.1016/j.cels.2017.01.010
- **[Tier 3, method definitions]** GIMME — Becker SA, Palsson BO (2008), *Context-Specific
  Metabolic Networks Are Consistent with Experiments*, PLoS Comput Biol 4(5):e1000082. DOI:
  10.1371/journal.pcbi.1000082. E-Flux — Colijn C et al. (2009), *Interpreting Expression Data
  with Metabolic Flux Models...*, PLoS Comput Biol 5(8):e1000489. DOI:
  10.1371/journal.pcbi.1000489. iMAT — Shlomi T et al. (2008), *Network-based prediction of
  human tissue-specific metabolism*, Nat Biotechnol 26(9):1003–1010, DOI: 10.1038/nbt.1487;
  Zur H, Ruppin E, Shlomi T (2010), *iMAT: an integrative metabolic analysis tool*,
  Bioinformatics 26(24):3140–3142, DOI: 10.1093/bioinformatics/btq602. pFBA — Lewis NE et al.
  (2010), *Omic data from evolved E. coli are consistent with computed optimal growth...*, Mol
  Syst Biol 6:390. DOI: 10.1038/msb.2010.47
- **[Tier 1]** *Colletotrichum* hemibiotrophy involves stage-specific effector/CAZyme
  secretion programs (real, but energetically opaque to this GAM category) — O'Connell RJ
  et al. (2012), *Lifestyle transitions in plant pathogenic Colletotrichum fungi deciphered by
  genome and transcriptome analyses*, Nat Genet 44(9):1060–1065. DOI: 10.1038/ng.2372
- **[Tier 2]** Appressorial penetration is a high-energy, turgor-driven process in related
  filamentous pathogens — Howard RJ et al. (1991), *Penetration of hard substrates by a fungus
  employing enormous turgor pressures*, PNAS 88(24):11281–11284, DOI: 10.1073/pnas.88.24.11281;
  de Jong JC et al. (1997), *Glycerol generates turgor in rice blast*, Nature 389:244–245,
  DOI: 10.1038/38418

## 4. Caveats

- **This is the central caveat, not a footnote.** The maximal method_range = 1.0 is a
  **structural artifact** of (a) the maintenance reaction lacking a GPR (no expression signal)
  and (b) the methods differing in whether they minimize flux. It should **not** be
  interpreted as biological variation in maintenance energetics.
- **Single-transcriptome limitation compounds it.** Only one static transcriptome (S1,
  pathogenic state, 3 technical reps) is applied to all media; expression does not vary by
  stage. Even if the maintenance reaction *had* a GPR, it could not report stage-specific
  energetics under this design.
- **Degenerate category size.** With `n_rxn = 1`, the "pathway" range reduces to a single
  reaction's ON/OFF disagreement, mechanically inflating it to the extremum — not comparable
  to multi-reaction pathways where range averages over many reactions.
- **The iMAT flux (~0.001) is numerically incidental**, not an estimated maintenance ATP
  demand; it reflects an unconstrained variable in the MILP, not a fitted NGAM value.
- **Fixable, if maintenance biology is of interest.** One could (i) fix ATPM to an
  empirically estimated lower bound so all methods carry the same maintenance flux, and/or
  (ii) exclude no-GPR reactions from method-concordance scoring. Either removes the artifact.
- **Alternative (weaker) reading:** iMAT's tolerance of low-cost cyclic/maintenance flux vs.
  the parsimony bias of the others is a *real* methodological difference — but it is a
  property of the solvers, not of *C. sublineola* infection energetics.

## 5. Confidence

**High** — that this specific discordance is a **modeling/method artifact**. The mechanism is
directly evidenced by the data (single no-GPR reaction; invariant binary split across all 18
conditions; iMAT-only nonzero flux at a fixed ~0.001) and is fully explained by well-
established GAM/NGAM and expression-integration conventions. Confidence that GAM discordance
reflects any *biological* infection signal is **Low** (essentially none, by construction).

```json
{
  "pathway": "GAM / maintenance (growth-associated & non-growth ATP maintenance)",
  "one_liner": "Maximal method disagreement (range=1.0) is a modeling artifact: the single no-GPR maintenance reaction gets no expression signal, so iMAT leaves it ON (~0.001) in all 18 conditions while parsimony-based pFBA/E-Flux/GIMME force it OFF — invariant to medium and stage, hence non-biological.",
  "confidence": "High (artifact); Low (biological signal)",
  "top_refs": [
    {"cite": "Thiele I, Palsson BO (2010) A protocol for generating a high-quality genome-scale metabolic reconstruction. Nat Protoc 5(1):93-121", "url": "https://doi.org/10.1038/nprot.2009.203"},
    {"cite": "Machado D, Herrgard M (2014) Systematic Evaluation of Methods for Integration of Transcriptomic Data into Constraint-Based Models of Metabolism. PLoS Comput Biol 10(4):e1003580", "url": "https://doi.org/10.1371/journal.pcbi.1003580"},
    {"cite": "Feist AM, Palsson BO (2010) The biomass objective function. Curr Opin Microbiol 13(3):344-349", "url": "https://doi.org/10.1016/j.mib.2010.03.003"},
    {"cite": "Orth JD, Thiele I, Palsson BO (2010) What is flux balance analysis? Nat Biotechnol 28(3):245-248", "url": "https://doi.org/10.1038/nbt.1614"},
    {"cite": "O'Connell RJ et al. (2012) Lifestyle transitions in plant pathogenic Colletotrichum fungi. Nat Genet 44(9):1060-1065", "url": "https://doi.org/10.1038/ng.2372"}
  ]
}
```
