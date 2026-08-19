# M3 Deep-Research Brief — FSP237 Infection-Stage Flux

## Organism & model
- **FSP237** = *Colletotrichum sublineola*, the causal agent of **sorghum anthracnose**
  (a hemibiotrophic fungal pathogen: an early **biotrophic** phase inside living host
  cells, switching to a destructive **necrotrophic** phase).
- **V10** genome-scale metabolic model (GEM), 1622 reactions.

## What the analysis did (context for interpretation)
We compared **vanilla FBA (pFBA)** against three **transcriptomics-constrained** flux
methods (**E-Flux, GIMME, iMAT**) across an **18-condition media panel** grouped by
anthracnose **infection stage** (pre-infection / biotrophic / necrotrophic / cocktail).

**Critical limitation (state in every interpretation):** there is only ONE transcriptome
(**S1**, the pathogenic state, 3 technical replicates), applied as a *static expression
prior* to every medium. Expression does NOT vary by stage — only the medium (hence flux)
does. So stage-differential signals are **flux/medium-driven**, not expression-driven.
No differential expression, no p-values on stages.

## Evidence tiers to use
- **Tier 1** — direct evidence in *Colletotrichum* spp. / anthracnose diseases.
- **Tier 2** — evidence in related hemibiotrophic / filamentous plant pathogens
  (*Magnaporthe oryzae*, *Fusarium*, *Zymoseptoria*, *Ustilago*, etc.).
- **Tier 3** — general fungal / biochemical inference (model organisms, textbook pathway
  logic) when no phytopathogen-specific evidence exists.

## Output contract (each pathway agent)
Return a self-contained markdown section with:
1. **Finding** — restate the quantitative result (numbers provided in the assignment).
2. **Mechanistic interpretation** — what the flux pattern plausibly means for infection.
3. **Evidence** — bulleted, each tagged [Tier 1/2/3], each with a real citation
   (author/year + DOI or URL you actually verified via search; do NOT fabricate).
4. **Caveats** — model/method artifacts, the single-transcriptome limitation, alternatives.
5. **Confidence** — High/Medium/Low with one-line justification.
End with a fenced ```json block: {"pathway","one_liner","confidence","top_refs":[{"cite","url"}]}.

Prioritize accuracy over completeness. If you cannot verify a citation, omit it and say so.
Do NOT invent DOIs or paper titles.
