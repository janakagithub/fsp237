# `fsp237_minimal_glucose/` — minimal-glucose growth test + ETC repair

The first cut at a glucose-viable FSP237 GSM, **predating** the
`fsp237_atp_safe_gsm` / biomass-extension / gap-fill pipeline. This is the
foundational notebook that:

1. Fetches the raw FSP237 model from KBase workspace `169876/251/1`
2. Adds all reactions from the curated *C. higginsianum* Excel set
3. Adds the minimum iMM904-essential reactions needed for biomass growth
   on glucose-minimal media at uptake = 5 mmol/gDW/h
4. Runs a final lower-MS-wins dedup pass

Plus a series of follow-up scripts that fixed ETC + futile cycle issues
discovered in the initial build.

## Files in this directory

### Build

| File | What it is |
|---|---|
| **`BuildMinimalFSP237Model.ipynb`** | The driving notebook. Produces `fsp237_minimal_glucose.json`. Run cells top-to-bottom. |
| `fsp237_initial_snapshot.json` | The pristine KBase fetch (workspace 169876/251/1) — input to the build. |
| **`fsp237_minimal_glucose.json`** | The final model after all the fixes below — biomass-viable on glucose, ETC works, futile cycles broken. |
| `fsp237_minimal_glucose_fluxes.tsv` | Per-reaction FBA flux dump at biomass = 0.226 (glc=−5). |
| `fsp237_minimal_glucose_flux_map.html` | Compartment-coloured Escher map of those fluxes (cytosol pink, mito green, ER blue, peroxisome orange, etc.). Open in a browser. |
| `fsp237_build_log.tsv` | Per-reaction add/skip/replace record from the build. |
| `fsp237_dedup_removed_log.tsv` | Reactions removed by the lower-MS-wins dedup pass and which canonical kept them. |
| `fsp237_dedup_pairs.xlsx` | Pretty Excel view of the dedup pairings. |
| `fsp237_excel_added_new.tsv` | Reactions added from the *C. higginsianum* Excel set. |
| `fsp237_imm904_essentials_kept.tsv` | Minimum iMM904-essential reactions retained. |

### Reaction de-duplication (first pass)

| Script | What it does |
|---|---|
| **`dedup_initial_build.py`** | Standalone CLI extraction of the dedup pass from `BuildMinimalFSP237Model.ipynb` (cell 15). Exact-stoichiometry hash (sign-sensitive, H+ ignored), lower-numeric-ModelSEED-ID-wins precedence, conservative GPR merge (adopt only if keeper has no GPR), bound widening, orphan-metabolite cleanup. Usage: `python dedup_initial_build.py <model.json> [<output.json>] [--log <log.tsv>]`. |

There is a **second, more sophisticated dedup script** for the final
gap-fill model — see `../simulations/gapfill_v1_v2/build_v3v4_dedup.py`.
The differences:

| | `dedup_initial_build.py` (this directory) | `build_v3v4_dedup.py` (later stage) |
|---|---|---|
| Stoichiometry hash | exact, sign-sensitive | both forward AND reverse hashes (collapses reverse-written duplicates) |
| GPR merge | conservative (only if keeper has none) | full union across all duplicates |
| Protected reactions | none | `bio1` (ATP-yield objective) |
| Compartment check | implicit via cpd IDs | explicit compartment-set match before merging |
| When to use | initial KBase+Excel+iMM904 merge, lots of same-direction duplicates | post-gap-fill cleanup, catches stoich-flipped pairs |

Use this script during/after the notebook build; use `build_v3v4_dedup.py`
on the gap-fill output.

### ETC / futile-cycle fixes (applied iteratively after the initial build)

| Script | What it fixes |
|---|---|
| **`add_complex_i.py`** | Adds a proton-pumping Complex I (`frxn13726_m0`): NADH·m + Q·m + 5 H⁺·m → NAD·m + QH₂·m + 4 H⁺·c |
| **`break_futile_cycles.py`** | Blocks 4 known ATP-loop futile cycle pairs (SAICAR alt-pathway, acyl-ACP bacterial cycle, ADP-sulfate cycle, NADH dehydrog. stoich-flip) by setting bounds to (0, 0) on the incorrect duplicates. |
| **`fix_etc_and_rerun.py`** | Locks frxn08975_c0 (External alt. NADH dehyd., inverted-stoich duplicate of rxn09524) and sets NGAM = 1.0. |
| **`force_etc_through_nadh.py`** | Eliminates 3 NADH-sink "cheats" (mito ethanol DH, malate DH reverse, citrate synthase reverse) so ETC has to be driven by real NADH oxidation. |
| `fsp237_6_direction_fixes_log.tsv` | Six initial-pass direction fixes from the notebook. |
| `fsp237_direction_fixes_log.tsv` / `_v2_log.tsv` | Later direction-audit fixes (irreversibilities). |
| `fsp237_er_vlcfa_route_log.tsv` | Notes on the ER VLCFA elongation route audit. |
| `fsp237_biomass_edits_proposed.md` / `.docx` | Working doc of biomass-coefficient adjustments proposed during this iteration (most were later carried into `fsp237_biomass_extension/`). |

## Workflow

```
KBase 169876/251/1 ──► BuildMinimalFSP237Model.ipynb ──► fsp237_minimal_glucose.json (v0)
                                                              ↓
                                                       add_complex_i.py
                                                              ↓
                                                       fix_etc_and_rerun.py
                                                              ↓
                                                       break_futile_cycles.py
                                                              ↓
                                                       force_etc_through_nadh.py
                                                              ↓
                                                       fsp237_minimal_glucose.json (final)
```

## Numbers at the end of this stage

- Biomass on glucose-minimal media (glc=−5, O₂=−10, NGAM=1.0): **0.226 1/h**
- ATP synthase (`rxn08173_m0`) flux: **0.997** — fully coupled OXPHOS
- O₂ uptake at biomass-max: −0.918 (73% increase over the initial-build's −0.532, reflecting real respiratory demand)
- NDI (mito NADH dehyd., non-pumping) flux: +0.787 — driving the chain
- Complex III flux: +0.811, Complex IV flux: +0.405

## What happens next

`fsp237_minimal_glucose.json` was the canonical FSP237 model before the
ATP-safe rebuild. The `../fsp237_atp_safe_gsm/` pipeline then built a
**cleaner** model from the curated CMM upward, using
`fsp237_minimal_glucose.json` only for media-bounds reference. The current
canonical model is `../simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version10_*.json`,
not this one — but this stage's output is still referenced by `extend_biomass.py`
as the source of GSM-style media bounds and as the source of additional
reactions/metabolites copied into the curated chain.
