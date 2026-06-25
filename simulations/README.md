# `simulations/` — multi-carbon-source simulations + gap-fill V1→V10

**Stage 4** of the FSP237 model-build pipeline. The GPR-updated model from
stage 3 is tested across **18 biologically motivated carbon-source
conditions** (each run aerobic + anaerobic = 36 simulations) representing
the three infection stages of *Colletotrichum sublineola* on sorghum
(pre-infection, biotrophic, necrotrophic), plus 3 stage-mimicking
cocktails. Failing conditions drive the gap-fill chain (V1 → V10).

## Top-level files

| File | What it is |
|---|---|
| **`run_simulation_panel.py`** | The simulation driver. 18 conditions × 2 O₂ states = 36 FBA runs. Reads a model JSON (configurable at the top), applies inorganic minimal media + condition-specific C-sources, optimises `bio_gsm`, dumps biomass + active fluxes + per-condition flux distribution. |
| **`simulation_results.tsv`** | Output — 36 rows with biomass, active flux count, C-uptake, biomass-per-C, solver status. The site Simulations tab pivots this to 18 condition rows with aer + ana side-by-side. |
| `condition_literature.tsv` | Per-condition literature support: 40 PMID-linked citations across 18 conditions, with one-line key findings. The site renders this as the **Literature** column on the Simulations tab. |
| **`RESULTS.md`** | Interpretation of the panel: which conditions grow, which fail, and the per-failure gap-fill recipe. |
| `per_condition/` | Per-condition non-zero flux dumps (≈25 files, one per condition × O₂ state where biomass > 0). |

## `gapfill_v1_v2/` — the V1 → V10 gap-fill chain

The 5 failing conditions in the initial panel (palmitate, oleate,
hexacosanoate, galacturonate, L-arabinose) drove a sequence of gap-fill
additions, each saved as its own model version. **All 10 versions are
preserved on disk** for reversibility and comparison.

### Build scripts (run in this order)

| Script | Builds | What it adds |
|---|---|---|
| **`build_v1_gapfill.py`** | V1 (`models/fsp237_gapfilled_Version1_noGenes.json`) | 28 new reactions: peroxisomal β-ox C8/C6/C4 closure + cofactor shuttles, Penttilä L-Ara pathway, Ashwell D-galU pathway. No GPRs. |
| **`build_v2_integrate_genes.py`** | V2 (`Version2_gapfill_genes_integrated.json`) | V1 + BLAST-derived FSP237 GPRs (19 reactions with confident hits) |
| **`build_v3v4_dedup.py`** | V3 / V4 | Exact-duplicate cleanup (33 dup groups → 34 rxns removed; `bio1` protected). GPRs union-merged across collapsed duplicates. Handles BOTH forward and reverse stoichiometry hashes (collapses reverse-written duplicates). See `../../fsp237_minimal_glucose/dedup_initial_build.py` for the simpler first-pass dedup used during the initial KBase + Excel + iMM904 merge. |
| **`build_v5v6_dirlock.py`** | V5 / V6 | β-ox + Ashwell direction-locked to degradation-only (34 reactions tightened). Prevents the optimizer from using these reactions in reverse for biosynthesis. |
| `build_v7v8_vlcfa.py` | V7 / V8 | First-pass VLCFA chain extension (had a cpd-id bug — superseded by V9/V10) |
| **`build_v9v10_vlcfa_cleanup.py`** | V9 / V10 | Complete VLCFA chain (C26→C24→C22→C20→C18, biosynthesis-consistent cpd IDs) + dropped 88 reactions in g0/n0/v0 compartments (zero flux across all 36 simulations) |

### Other tools

| File | What it is |
|---|---|
| **`find_candidates.py`** | Phase 2+3 of the gap-fill workflow: pulls candidate proteins from `C_higgensium.gbff` by EC + product-name scoring, BLASTs against the FSP237 DB, picks the best per gap-fill reaction. Output → `blast/v1_rxn_to_gene_mapping.tsv`. Used by `build_v2_integrate_genes.py`. |
| **`test_gapfilled_model.py`** | Re-run the 18-condition panel against any model. Usage: `python test_gapfilled_model.py models/<model>.json v<N>` → results under `reports/v<N>_simulation_results.tsv`. |
| **`save_v6_to_kbase.py`** | Push any V* model to KBase workspace 28277 as a `KBaseFBA.FBAModel` object with proper compartments + biomasses + GPR feature_refs that resolve in the *C. sublineola* JGI genome 169876/166/3. |
| `faoxidation_escher_map.tiff` | Screenshot of the peroxisomal FA-oxidation area on the existing Escher map (used during the VLCFA gap-fill diagnosis). |

### Subdirectories

- **`models/`** — 10 versioned model JSONs (V1 through V10). Each ~1MB. Load with `cobra.io.load_json_model(path)`.
- **`candidates/`** — C. higginsianum candidate FASTAs (`v1_query_candidates.faa`) + provenance TSV.
- **`blast/`** — Raw BLAST output + final rxn→FSP237_gene mapping table (`v1_rxn_to_gene_mapping.tsv`).
- **`reports/`** — Per-version panel results, dedup/dirlock/VLCFA logs, biological narrative docs:
    - **`SUMMARY.md`** — master narrative of every V1→V10 decision
    - **`PATHWAY_DIAGRAMS.md`** + `.docx` — full diagrams of β-ox + Penttilä + Ashwell pathways with equations
    - **`OLEATE_PATHWAY.md`** + `.docx` — focused oleate degradation pathway diagram
    - **`OLEATE_LITERATURE_SUPPORT.md`** + `.docx` — literature justification for oleate as a pre-infection carbon source
    - `v{1..10}_simulation_results.tsv` — 36-row panel results per version
    - `v{1..10}_per_condition/` — per-condition flux dumps per version
    - `v3_dedup_log.tsv`, `v4_dedup_log.tsv` — dedup decision audit
    - `v5_direction_locks.tsv` — per-rxn before/after bounds for the dirlock pass
    - `v7_vlcfa_added.tsv`, `v9_change_log.tsv` — VLCFA gap-fill audit
    - `v1_added_reactions.tsv` — V1 gap-fill reactions with status/EC/equation/reason
    - `v2_gpr_assignments.tsv` — V2 BLAST gene assignment decisions per reaction

## Workflow

```
fsp237_atp_safe_gsm_gpr_updated.json (from ../gpr-update/)
        ↓
        ├─► run_simulation_panel.py  (panel results on the source model)
        │
        └─► build_v1_gapfill.py       → V1 (no genes)
                ↓
            build_v2_integrate_genes.py  → V2 (with genes)  [also: find_candidates.py first]
                ↓
            build_v3v4_dedup.py       → V3 / V4 (dedup applied to V1 / V2)
                ↓
            build_v5v6_dirlock.py     → V5 / V6 (β-ox + Ashwell dirlocked)
                ↓
            build_v9v10_vlcfa_cleanup.py  → V9 / V10 (complete VLCFA + compartment prune)
                ↓
            test_gapfilled_model.py models/Version10... v10  → reports/v10_*
                ↓
            save_v6_to_kbase.py       → KBase 28277/287/1
```

The current canonical / published model is **V10**:
`models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json`.

## Headline results (V10, aerobic)

| Stage | Best condition | Biomass (1/h) |
|---|---|---:|
| Necrotrophic (cocktail) | suc + maltose + galU + xyl + ara + 20 AA | **0.446** |
| Necrotrophic (single) | maltose (starch proxy) | 0.364 |
| Pre-infection (single) | hexacosanoate (C26) | 0.206 |
| Biotrophic (cocktail) | glc + suc + Glu + Gln + malate | 0.063 |
| Pre-infection (cocktail) | trehalose + glycerol + oleate | 0.069 |
| Biotrophic (single) | L-Glutamate / L-Glutamine | 0.053 |

All 18 aerobic conditions grow in V10. Anaerobic FA / GABA / L-Ara
conditions correctly fail (no ETC means no NADH/FADH₂ reoxidation —
biological, not bug).

## Dependencies

- Python 3.10+ with `cobra`, `pandas`, `subprocess` (BLAST)
- `blastp` on PATH (for `find_candidates.py`)
- `cobrakbase` (for `save_v6_to_kbase.py` only)
- KBase token at `~/.kbase/token` (for the save script)
