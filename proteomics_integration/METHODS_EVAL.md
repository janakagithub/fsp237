# Proteomics × FBA — methods evaluation + integration plan (Phases 0–2)

**Scope.** How to use the FSP237 proteomics (`proteomics-fsp237/`) inside the V10
GEM. Proteomics-first; metabolomics is a deliberate later add. Evaluates the five
requested methods (GECKO/ecGEM, sMOMENT, IOMA, TFA/TMFA, pyTFA) against the data we
actually have, then lays out a concrete Phase 0–2 build.

> Web verification was unavailable in the authoring session (org policy blocked
> search). DOIs supplied by the requester are trusted as given; DOIs added below are
> from memory and should be confirmed before manuscript use. Tooling API specifics
> flagged **[verify]** should be checked against current docs at build time.

---

## 1. The data that gates method choice

Two datasets under `proteomics-fsp237/`:

- **`Cs_MPLEx_proteins_mycoCosum/` — the FBA-critical set.** 4,981 Cs proteins,
  native **`gene_NNNN`** IDs, across **PDA / half-PDA / one-tenth-PDA, 5 reps each**
  (15 samples), TMT **log2 abundance** (~19–36). `protein_wide_filtered_imputed.csv`
  has **0 NaN** (analysis-ready); the non-imputed file is ~18% missing. Three limma
  DE contrasts (half-vs-PDA = 150 sig, onetenth-vs-PDA = 1,204, onetenth-vs-half =
  1,237 at adjP<0.05). **Coverage of V10:** 832/1,274 model genes (65%), touching
  **882/1,045 GPR reactions (84%)**.
- **`global_proteomics/BRAVE_*` — in-planta dual proteome.** 7,764 proteins, UniProt
  IDs, sorghum `SORBI_*` host + Cs, genotypes BT623 (susceptible) / SC112 (resistant)
  × control/infected. Host–pathogen biology, **not** a media condition — a separate
  track, not a PDB-media FBA input.
- **`annota_final3_new.xlsx` — the bridge.** Per `gene_NNNN`: `uniprot_chig_id`,
  **`ko_id`**, `kegg_hit`, `pfam_hits`, C. higginsianum description, "identified in
  proteomics?", "unique to one-tenth PDA". Links `gene_ ↔ UniProt ↔ KEGG/EC`.

Four properties decide feasibility:

| Property | Value here | Consequence |
|---|---|---|
| Quantification | **Relative** (TMT log2 ratios) | Per-enzyme absolute capacity bounds `v≤kcat·[E]` need a relative→absolute assumption; **pool-level** constraints do not. |
| Metabolomics | **None (for now)** | Hard blocker for IOMA; removes the tightening role of TFA. |
| kcat set for Cs | **None** experimental | Must parameterize (BRENDA/SABIO hierarchy + ML predictors). |
| IDs / condition match | `gene_NNNN`, 84% GPR coverage; **proteome conditions = the 3 simulated PDB media** | Direct GPR mapping and genuine cross-condition integration — the project's main asset. |

---

## 2. Method-by-method verdict

| Method | Consumes proteomics? | Needs metabolomics? | Needs kcat? | Needs Km? | Needs *absolute* protein? | Fit now |
|---|---|---|---|---|---|---|
| **GECKO / ecGEM** | ✅ core | no | yes | **no** | ideally; workable via pool | **★ Primary** |
| **sMOMENT** | pool-level | no | yes | **no** | no (single pool) | **★ Fast baseline** |
| **IOMA** | ✅ (paired) | **yes — blocker** | yes | **yes** | yes | ⛔ Deferred |
| **TFA / TMFA / pyTFA** | ❌ (thermodynamic) | optional | no | no | n/a | ◐ Later / complementary |

### GECKO / ecGEM — right center of gravity
Purpose-built for "proteomics → enzyme capacity." Degrades gracefully: trusted
enzymes get individual `v ≤ kcat·[E]` bounds; the rest draw from a shared
**total-protein pool** (`Σ MWₑ·[E]ₑ ≤ P·f·σ`) needing only a scalar total-protein
content. With **relative** TMT this means: build one kcat-parameterized ecModel
skeleton, then per condition **allocate the enzyme budget by the measured proteome**
and test whether that reshapes flux and reproduces the biomass dilution response
better than uniform pFBA. Km is **not** used (capacity constraint assumes saturation).
Lift = kcat parameterization + a defensible relative→pool mapping.
**[verify]** Python `geckopy` API vs the more mature MATLAB GECKO3 toolbox.

### sMOMENT (+ AutoPACMEN) — tractable first deliverable
Enzymes lumped into a single pool; AutoPACMEN automates kcat/MW retrieval; needs **no
absolute per-protein data**. Trade-off: collapses 4,981 per-condition measurements to
one pool bound, under-using the richest signal. Ideal **week-1 baseline** that
de-risks the kcat pipeline before full GECKO; machinery transfers directly. Km not used.

### IOMA — genuinely blocked without metabolomics
Fits Michaelis–Menten rate laws in a QP: `v = f([E], [S], kcat, Km)`. Needs enzyme
abundance **and** metabolite concentrations **and** per-reaction Km+kcat. Missing the
metabolomics half means the rate equations can't be posed — not a soft limitation.
Natural payoff **when metabolomics arrives** (uniquely consumes both layers jointly).
No maintained COBRApy implementation → from-scratch build.

### TFA / TMFA / pyTFA — valuable but not a proteomics method
Thermodynamic constraints (ΔrG directionality, loopless feasibility) come from ΔfG
(group/component contribution), **not** proteomics — so it doesn't answer "how to use
the proteomics." Complementary, and **ModelSEED-native** (matches our `cpd*` IDs;
pyTFA ships a group-contribution thermo DB — **[verify]** seed-ID coverage of our
metabolite set). Grows more powerful with measured metabolite ranges. Best role:
layer onto the finished ecModel later (= **ETFL**, enzyme+thermo).

---

## 3. How kcat and Km are obtained (the parameterization question)

**First, the clarifying fact:** for the proteomics-first enzyme-constrained methods we
would actually run now — **GECKO, sMOMENT, MOMENT — only kcat is needed; Km is not.**
They are *capacity* constraints (`v ≤ kcat·[E]`) that assume enzymes operate at their
turnover number (saturated), so no saturation/affinity term enters. **Km only becomes
relevant for IOMA** (and any true kinetic/rate-law model), which is deferred to the
metabolomics phase. Thermodynamic TFA/pyTFA uses **neither** kcat nor Km (it uses ΔfG).

### 3.1 kcat — four sources, used in a cascade

**Per reaction we need one kcat (the catalyzing enzyme's turnover, 1/s → ×3600 = 1/h)
and each enzyme's MW (g/mmol, summed from its sequence).** For C. sublineola there is
no experimental kinetics, so:

1. **Curated-DB hierarchy (BRENDA / SABIO-RK).** Match on EC number with graceful
   fallback: `(EC, organism, substrate)` → `(EC, substrate)` → `(EC, any organism)` →
   wildcard EC. Take the **maximum** across matches (measured Vmax reflects best-case;
   GECKO uses kcat as an upper capacity). **We already have EC for ~794/1,622 reactions
   (49%)** from the earlier KEGG mapping (`rna_seq_integration/.../data/rxn_kegg_map.tsv`)
   plus `ko_id`→EC via the annotation xlsx — a strong starting index.
2. **ML prediction — the primary fill for a data-poor fungus.**
   - **DLKcat** (Li 2022): predicts kcat from **enzyme amino-acid sequence + substrate
     SMILES**. We have both — sequences in `.../data/fsp237_proteome.faa` (14,857
     proteins) and substrate SMILES from the local ModelSEEDDatabase (`cpd*`→SMILES/
     InChI). Per reaction: GPR→gene→sequence, substrate = principal non-cofactor
     reactant; predict per (enzyme, substrate) pair; aggregate (max) per reaction.
     GECKO 3.0 integrates DLKcat directly.
   - **TurNuP** (Kroll 2023): kcat from enzyme sequence + a reaction fingerprint;
     avoids explicit enzyme–substrate pairing, useful where the principal substrate is
     ambiguous. Good cross-check on DLKcat.
3. **kapp / kmax from our own data (refinement/sanity, not primary).** Apparent
   catalytic rate `kapp,i = v_i / [E]_i`; `kmax,i = max` over conditions (Davidi 2016;
   Heckmann 2020). Uses the 3 dilutions as the condition set. **Requires absolute
   enzyme**, so with relative TMT this is a *consistency check* on DB/ML kcat, not a
   source — unless we adopt the pseudo-absolute conversion (§4).
4. **GECKO automatic kcat tuning.** After assembly, if predicted µmax < observed,
   raise the kcat of the enzyme with highest flux control (it was likely the
   underestimate). GECKO3 `sensitivityTuning`; a small manual curation loop otherwise.

**Assignment cascade we will use:** experimental (BRENDA/SABIO exact-ish) → DLKcat/
TurNuP prediction → EC-family median → global default; then GECKO tuning on the
growth-limiting subset. Every reaction's kcat provenance is logged (a `kcat_source`
column) so the manuscript can report the fraction experimental vs predicted.

### 3.2 Km — only for the deferred IOMA/metabolomics phase

When metabolomics lands: Km from **BRENDA / SABIO-RK** (same EC hierarchy as kcat), or
**ML Km predictors** (e.g. Kroll et al. Km models) from sequence + substrate, paired
with the measured metabolite concentration `[S]`. Not needed for anything in Phases 0–3.

### 3.3 MW and units (needed regardless)

- **MW** per enzyme from its sequence (Biopython `ProtParam`, sum of residue masses),
  g/mmol.
- **Enzyme-usage stoichiometry (GECKO):** a reaction with flux `v` (mmol/gDW/h) draws
  its enzyme at `v/kcat` (mmol enzyme/gDW), i.e. coefficient `1/kcat` on an
  enzyme-usage pseudo-reaction; the pool constraint sums `MWₑ·[E]ₑ ≤ P·f·σ`
  (g protein/gDW). `P` = total protein content; `f` = mass fraction that is
  enzyme/modeled; `σ` = average saturation. kcat converted 1/s → 1/h (×3600).

---

## 4. The one decision to make before any enzyme-constrained build

**How to treat relative TMT abundance.** Two options:

- **(A) Total-protein pool + relative allocation (recommended, assumption-light).**
  Keep abundances relative; use a single pool bound and let the per-condition proteome
  set *relative* enzyme availability. Tests the real question — does reallocating a
  fixed enzyme budget per the measured proteome reshape flux across the dilution
  series — without inventing absolute numbers.
- **(B) Pseudo-absolute conversion (enables per-enzyme bounds; heavier assumptions).**
  Convert TMT→absolute via total-protein anchoring / iBAQ-style rescaling. Unlocks
  individual `v≤kcat·[E]` bounds and kapp estimation, but TMT reporter intensities
  support this weakly. Revisit if label-free/iBAQ intensities become available.

Default to **(A)**; document the choice in the payload `meta`.

---

## 5. Integration plan — Phases 0–2 (near term)

Reuses the frozen RNA-seq machinery: `rna_seq_integration/src/gpr_expression.py`
(`aggregate`, AND=min / OR=max), `stage2_gimme_imat_eflux.py`, `stage1_overlay.py`,
`stage5_pathway_analysis.py`, and the `concordance()` metric. New module lives in
`proteomics_integration/` mirroring the RNA-seq layout. **No model, media, or biomass
change; CPLEX/pFBA reused as in the PDB-media panel.**

Condition mapping is 1:1:

| Proteomics condition | FBA media key (`run_simulation_panel.py`) |
|---|---|
| `fullstrengthPDA_{1..5}` | `19_pdb_baseline` |
| `halfPDA_{1..5}` | `20_pdb_half` |
| `onetenthPDA_{1..5}` | `21_pdb_onetenth` |

### Phase 0 — Ingest, QC, harmonize
- Load `protein_wide_filtered_imputed.csv`; map columns→conditions; collapse 5 reps →
  per-condition mean (retain per-rep for bootstrap).
- Normalization QC (column medians/distribution); confirm log2 scale.
- Build `gene_ → {uniprot, ko, EC, pfam, desc, detected, unique_onetenth}` from
  `annota_final3_new.xlsx`.
- Coverage report vs V10 GPRs (per-condition detection; the 832/882 baseline).
- **Outputs:** `outputs/proteome_condition_means.tsv`, `outputs/gene_annotation_map.tsv`,
  `outputs/coverage_report.md`.

### Phase 1 — Condition-matched proteome-constrained flux
- Per-reaction per-condition enzyme score via `aggregate` (AND=min, OR=max).
- Run **E-Flux, GIMME, iMAT** per condition (PDA/half/onetenth) on the matched PDB
  media, plus the **pFBA** baseline. Parameterize `stage2_gimme_imat_eflux.py` to take
  the proteome matrix instead of the single S1 vector.
- **Key novelty vs the S1 work:** expression now *varies by condition* — genuine
  cross-condition integration, not one fixed prior.
- **Outputs:** `flux_results/wide_flux_{pfba,eflux,gimme,imat}.tsv`,
  `flux_results/unified_flux_matrix.tsv` (mirror the infection-stage layout).

### Phase 2 — Concordance + dilution-response validation
- Per-condition concordance of each method's flux vs the proteome (reuse
  `concordance()`: Spearman |flux|~expr, hi-recall, precision, MCC/Jaccard, composite).
- **Dilution-response test (the payoff):** does predicted flux change PDA→onetenth
  track the measured protein **logFC direction** (the 3 DE tables as ground truth)?
  Directional-agreement scored per reaction and per pathway (23-bucket + KEGG).
- **Proteome vs S1-transcriptome cross-check:** does protein-constrained flux agree
  with the earlier transcript-constrained flux?
- **Outputs:** `method_comparison/concordance_*.tsv`,
  `statistics/dilution_response.tsv`, `RESULTS.md`.

### Deliverable framing
A self-contained "proteomics-constrained flux" milestone paralleling
`infection_stage_flux_analysis`, delivering the condition-matched validation that is
the whole point of this dataset — and the input that motivates the enzyme-constrained
track below.

---

## 6. Roadmap beyond Phase 2 (enzyme-constrained + later metabolomics)

- **Phase 3 — kcat parameterization.** EC index (from KEGG map + `ko_id`) → BRENDA/
  SABIO hierarchy; DLKcat/TurNuP on the proteome FASTA for gaps; MW from sequence;
  provenance-logged kcat table.
- **Phase 4 — sMOMENT/AutoPACMEN pool ecModel.** Baseline enzyme-constrained model;
  validates the kcat pipeline; per-condition pool bound.
- **Phase 5 — GECKO ecModel with per-condition proteome allocation.** Benchmark vs
  pFBA and the Phase 1 E-Flux/GIMME/iMAT results across the 3 dilutions.
- **Phase 6 (needs metabolomics) — IOMA and pyTFA→ETFL.** IOMA as the joint
  proteomics+metabolomics kinetic endpoint; pyTFA thermodynamic layer, then ETFL
  (enzyme+thermo) on the finished ecModel.

---

## 7. Caveats to carry into the manuscript

- PDB-media biomass response is trivially glucose-scaled; the scientific question is
  whether the **proteome constraint reshapes flux beyond that scaling** — state it
  explicitly and center Phase 2 on it.
- TMT log2 is **relative**; every absolute-flavored step (per-enzyme bounds, kapp)
  inherits the §4 assumption — label it.
- kcat for a fungus is mostly **predicted, not measured**; report the experimental vs
  predicted fraction and DLKcat/TurNuP agreement.
- Use the **imputed** matrix for modeling; note imputation of the one-tenth-PDA gaps.
- BRAVE (in-planta) is a **separate track** (host + pathogen; UniProt) — do not feed
  `SORBI_*` host proteins into the fungal model.

---

## 8. References

Requester-supplied (trusted as given):
- **GECKO / ecGEM** — Sánchez et al. 2017, *Mol Syst Biol* — DOI 10.15252/msb.20167411
- **IOMA** — Yizhak et al. 2010, *Bioinformatics* — DOI 10.1093/bioinformatics/btq183
- **TMFA** — Henry, Broadbelt & Hatzimanikatis 2007, *Biophys J* — DOI 10.1529/biophysj.106.093138
- **pyTFA** — Salvy et al. 2019, *Bioinformatics* — DOI 10.1093/bioinformatics/bty499
- **sMOMENT / AutoPACMEN** — Bekiaris & Klamt 2020, *BMC Bioinformatics* — DOI 10.1186/s12859-019-3329-9

Added (confirm DOIs before manuscript):
- GECKO 3.0 protocol — Domenzain et al. 2022, *Nat Commun* 13:3766 — DOI 10.1038/s41467-022-31421-1
- MOMENT — Adadi et al. 2012, *PLoS Comput Biol* — DOI 10.1371/journal.pcbi.1002575
- DLKcat — Li et al. 2022, *Nat Catal* — DOI 10.1038/s41929-022-00798-z
- TurNuP — Kroll et al. 2023, *Nat Commun* 14:4139 — DOI 10.1038/s41467-023-39840-4
- kapp/kmax — Davidi et al. 2016, *PNAS* — DOI 10.1073/pnas.1514240113
- ML kcat — Heckmann et al. 2020, *Nat Commun* — DOI 10.1038/s41467-020-19070-8
- ETFL — Salvy & Hatzimanikatis 2020, *Nat Commun* — DOI 10.1038/s41467-019-13818-7
- E-Flux — Colijn et al. 2009, *PLoS Comput Biol* — DOI 10.1371/journal.pcbi.1000489
- GIMME — Becker & Palsson 2008, *PLoS Comput Biol* — DOI 10.1371/journal.pcbi.1000082
- iMAT — Shlomi et al. 2008, *Nat Biotechnol* — DOI 10.1038/nbt.1487
