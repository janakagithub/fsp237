# Cell-wall polysaccharide metabolism

*(fungal cell-wall biosynthesis/remodeling — chitin, β-glucan, mannan, GPI-anchored wall proteins)*
Organism: **FSP237 = *Colletotrichum sublineola*** (sorghum anthracnose; hemibiotroph). Model: V10 GEM (1622 reactions).

## 1. Finding

Stage-differential flux (active-rate) for cell-wall polysaccharide metabolism:

| Infection stage (media group) | Active-rate |
|---|---|
| Pre-infection | 0.4444 |
| **Biotrophic** | **0.4778** (highest) |
| Necrotrophic | 0.4333 |

- **Stage range = 0.0445** — modest; cell-wall flux is slightly elevated in biotrophic-stage media and lowest in necrotrophic-stage media.
- **Method discordance:** mean_method_range = 0.1976 (moderate) across pFBA vs. E-Flux / GIMME / iMAT.

The ordering biotrophic > pre-infection > necrotrophic is directionally consistent with the known biology of the hemibiotrophic switch, but the absolute stage spread is small.

## 2. Mechanistic interpretation

The fungal cell wall — an inner load-bearing **chitin + β-1,3/β-1,6-glucan** core, an outer **mannan / α-1,3-glucan** and **GPI-anchored glycoprotein** layer — is continuously rebuilt as *Colletotrichum* transits spore → germ tube → melanized appressorium → **biotrophic infection vesicle/primary hyphae** → **necrotrophic secondary hyphae**. Each morphotype requires de novo polysaccharide synthesis (chitin synthases, the FKS β-1,3-glucan synthase, mannosyltransferases) coordinated with lytic/remodeling enzymes (chitinases, glucanases, chitin deacetylases) under cell-wall-integrity (CWI) MAPK control.

A modestly **higher biotrophic flux** is the biologically expected direction for two reasons:

1. **Immune-evasion "masking."** The biotrophic infection vesicle is enveloped by host plasma membrane and must avoid triggering chitin-triggered immunity (CTI). Phytopathogens actively shield the immunogenic chitin core: converting surface chitin to less-recognized **chitosan** via **chitin deacetylases**, and/or coating the wall with **α-1,3-glucan** that physically blinds host chitin receptors. Sustaining these outer-layer polymers during biotrophy imposes an anabolic (precursor) demand on cell-wall polysaccharide reactions — plausibly the small biotrophic bump seen here.
2. **Precursor availability, not induced expression.** Because only the medium (not the transcriptome) varies here, the biotrophic peak most parsimoniously reflects **medium-driven precursor supply** — flux carried by UDP-GlcNAc (chitin), UDP-glucose (β-glucan), and GDP-mannose (mannan) pools that the biotrophic-stage media happen to feed better than the nutrient-restricted necrotrophic set. It is *not* evidence of stage-induced wall-gene up-regulation.

In necrotrophy the fungus proliferates in dead/dying tissue and shifts investment toward secreted CWDEs and rapid hyphal biomass; a slightly lower normalized cell-wall-polysaccharide active-rate is consistent with that reallocation but, again, is medium-driven here.

**Virulence / fungicide context:** chitin synthases, the FKS β-1,3-glucan synthase (the echinocandin target), and CWI signaling are established antifungal/fungicide nodes. In *Colletotrichum* specifically, β-1,3-glucan synthesis (Rho4), chitin synthase CfCHS1, and CWI-MAPK signaling are each required for full virulence — making this pathway a credible target-context readout even when the stage signal is weak.

## 3. Evidence

**Tier 1 — *Colletotrichum* spp. / anthracnose**
- β-1,3-glucan synthesis and CWI are required for virulence in the **maize anthracnose** fungus *C. graminicola* — the closest well-studied relative of grass-infecting *C. sublineola*: the Rho4 GTPase governs β-1,3-glucan synthesis, wall integrity and full virulence. Corrêa dos Santos et al., 2022, *J. Fungi* 8:997. [Tier 1] https://doi.org/10.3390/jof8100997
- Chitin synthase **CfCHS1** mediates cell-wall integrity, stress tolerance and pathogenicity in *C. fructicola*. Liu et al., 2023, *J. Fungi* 9:643. [Tier 1] https://doi.org/10.3390/jof9060643
- The **CWI MAP-kinase pathway** is required for development, pathogenicity and stress adaptation in the anthracnose fungus *C. scovillei*. 2023, *Mycobiology* 51:1. [Tier 1] https://doi.org/10.1080/12298093.2023.2220171
- Transcription factor **Con7** regulates CWI, appressorium/hyphopodium formation and pathogenicity in *C. graminicola* and *C. siamense*. 2024, *J. Fungi* 10:495. [Tier 1] https://doi.org/10.3390/jof10070495
- DHN-melanin genes **ChPks/ChThr1** are needed for wall integrity and pathogenicity in the hemibiotroph *C. higginsianum* (links appressorial wall/melanin layer to penetration). 2023, *Int. J. Mol. Sci.* 24:15890. [Tier 1] https://doi.org/10.3390/ijms242115890
- Calcineurin-responsive TF **CgCrzA** is required for CWI and infection-related morphogenesis in *C. gloeosporioides*. 2020, *Plant Pathol. J.* 36. [Tier 1] https://doi.org/10.5423/PPJ.OA.04.2020.0071

**Tier 2 — related phytopathogenic / filamentous fungi**
- **Surface α-1,3-glucan enables fungal "stealth" infection** by masking cell-wall components from plant innate immunity — foundational chitin/β-glucan masking evidence in the hemibiotroph *Magnaporthe (Pyricularia) oryzae* and other pathogens. Fujikawa et al., 2012, *PLoS Pathog.* 8:e1002882. [Tier 2] https://doi.org/10.1371/journal.ppat.1002882
- Chitin deacetylase **PoCda7** contributes to pathogenicity of *Pyricularia oryzae* (chitin→chitosan conversion tied to virulence). 2021, *Microbiol. Res.* 251:126749. [Tier 2] https://doi.org/10.1016/j.micres.2021.126749
- Chitin deacetylation / chitosan generation during vegetative growth in *M. oryzae*. Geoghegan & Gurr, 2017, *Cell. Microbiol.* 19:e12743. [Tier 2] https://doi.org/10.1111/cmi.12743
- **Inhibiting chitin deacetylases attenuates plant fungal diseases** — validates wall-remodeling enzymes as fungicide targets. 2023, *Nat. Commun.* 14:3853. [Tier 2] https://doi.org/10.1038/s41467-023-39562-7
- GPI-anchored wall protein **FocECM33** regulates growth and virulence in *Fusarium oxysporum* f. sp. *cubense* TR4 (GPI-wall-protein layer in a filamentous phytopathogen). 2022, *Fungal Biol.* 126. [Tier 2] https://doi.org/10.1016/j.funbio.2021.12.005

**Tier 3 — general fungal cell-wall biology / drug-target biochemistry**
- Cryo-EM structure & mechanism of the fungal **β-1,3-glucan synthase FKS1** (echinocandin target). Hu et al., 2023, *Nature* 616. [Tier 3] https://doi.org/10.1038/s41586-023-05856-5
- Independent structure of a fungal 1,3-β-glucan synthase. 2023, *Sci. Adv.* 9:eadh7820. [Tier 3] https://doi.org/10.1126/sciadv.adh7820
- The fungal cell wall as a target for new antifungal therapies (chitin, β-glucan, mannan, GPI layers). Hopke et al., 2019, *Biotechnol. Adv.* 37. [Tier 3] https://doi.org/10.1016/j.biotechadv.2019.02.008
- "Chitins and chitosans — a tale of discovery and disguise, of attachment and attainment" — review of wall-carbohydrate masking/immune disguise. 2024, *Curr. Opin. Plant Biol.* 82:102661. [Tier 3] https://doi.org/10.1016/j.pbi.2024.102661

*(All citations retrieved and confirmed via NCBI PubMed / DOI; titles, journals and DOIs are as indexed. Note: PMID 39536646 was indexed without named first author in the summary record, so it is cited by title.)*

## 4. Caveats

- **Single-transcriptome limitation (decisive):** one static expression prior (S1 pathogenic state, 3 technical reps) is applied to every medium. Expression does not vary by stage, so the biotrophic peak is **medium/flux-driven, not evidence of stage-induced wall-gene expression**. Any "biotrophic masking" narrative is a mechanistic *plausibility* consistent with the direction of the flux, not a measured induction.
- **Small effect size:** stage range 0.0445 on an active-rate near ~0.45 is a ~10% relative spread — modest and easily within model/medium noise. Do not over-interpret the biotrophic > necrotrophic ordering.
- **Moderate method discordance (0.1976):** pFBA vs. E-Flux/GIMME/iMAT disagree moderately; the constrained methods weight the same static transcriptome differently, so the absolute active-rate is method-sensitive.
- **GEM abstraction:** cell-wall "polysaccharide metabolism" here is a reaction set with lumped biomass/precursor demands; it cannot distinguish chitin vs. β-glucan vs. mannan sub-fluxes, nor resolve outer-layer masking (α-1,3-glucan, chitin→chitosan) unless those reactions are explicitly encoded. Remodeling/lytic enzymes and CWI signaling are regulatory, largely outside stoichiometric flux.
- **Alternative explanation:** the biotrophic bump may simply track UDP-GlcNAc / UDP-glucose / GDP-mannose precursor availability in the biotrophic media rather than any wall-specific program.

## 5. Confidence

**Medium-Low.** The mechanistic story (biotrophic-phase chitin/β-glucan masking; CWI/chitin-synthase/FKS as virulence and fungicide nodes) is strongly supported by Tier-1 *Colletotrichum* and Tier-2 phytopathogen literature, but the *quantitative* stage signal is small (range 0.0445), medium-driven, and confounded by the single static transcriptome plus moderate method discordance — so the flux result corroborates known biology only weakly and cannot itself demonstrate stage-specific wall regulation.

```json
{
  "pathway": "Cell-wall polysaccharide metabolism (chitin / β-glucan / mannan / GPI wall proteins)",
  "one_liner": "Cell-wall biosynthetic flux is modestly highest in biotrophic-stage media (0.4778 vs 0.4444 pre / 0.4333 necro; range 0.0445), directionally consistent with biotrophic chitin/β-glucan masking for immune evasion, but the small, medium-driven signal (single static transcriptome, method range 0.1976) makes it corroborative rather than demonstrative.",
  "confidence": "Medium-Low",
  "top_refs": [
    {"cite": "Corrêa dos Santos et al., 2022, J. Fungi 8:997 — Rho4 β-1,3-glucan synthesis, cell-wall integrity and virulence in Colletotrichum graminicola (maize anthracnose)", "url": "https://doi.org/10.3390/jof8100997"},
    {"cite": "Fujikawa et al., 2012, PLoS Pathog. 8:e1002882 — Surface α-1,3-glucan enables fungal stealth infection by masking the wall from plant immunity", "url": "https://doi.org/10.1371/journal.ppat.1002882"},
    {"cite": "Liu et al., 2023, J. Fungi 9:643 — Chitin synthase CfCHS1, cell-wall integrity and pathogenicity in Colletotrichum fructicola", "url": "https://doi.org/10.3390/jof9060643"},
    {"cite": "2023, Nat. Commun. 14:3853 — Inhibition of chitin deacetylases to attenuate plant fungal diseases", "url": "https://doi.org/10.1038/s41467-023-39562-7"},
    {"cite": "Hu et al., 2023, Nature 616 — Structural and mechanistic insights into fungal β-1,3-glucan synthase FKS1 (echinocandin target)", "url": "https://doi.org/10.1038/s41586-023-05856-5"}
  ]
}
```
