# `fsp237_atp_safe_gsm/` — ATP-safe expansion of the FSP237 CMM into a GSM

This is **stage 1** of the FSP237 model-build pipeline. Goal: take the
curated central-carbon model (CMM) and iteratively expand it into a
genome-scale model (GSM) **without breaking the ATP yields** (30 ATP/mmol
glucose aerobic, 2 anaerobic).

The CMM gives exactly 30/2 by construction. The GSM contains many extra
reactions, some of which form futile ATP cycles. Naively unioning the two
inflates the yield to nonsense. This pipeline picks reactions one at a time
and rejects any that would push ATP yield above the curated baseline.

## Files in this directory

| File | What it is |
|---|---|
| **`BuildFSP237_ATPSafe_GSM.ipynb`** | The driving Jupyter notebook. Run cells top-to-bottom. |
| **`atp_safe_expand.py`** | Helper module: ATP-yield testing, candidate expansion in bisected batches, shared-reaction swap search, blacklist re-introduction. The notebook calls into this. |
| `REPORT.md` | Narrative summary of the build. Read this first if you don't want to re-run the notebook. |
| `cmm_169876_170_5.json` | Cached KBase fetch of the curated CMM (workspace `169876/170/5`). The notebook re-fetches if this is missing. |
| **`fsp237_atp_safe_gsm.json`** | The output of this stage — the ATP-safe GSM. Downstream stages (biomass extension, GPR update, simulations) all start from this file. |
| `safe_added.tsv` | Every reaction added safely during the expansion. |
| `blacklist.tsv` | Every reaction rejected because it would have pushed ATP yield above 30/2. |
| `shared_rxns_swapped.tsv` | Reactions that needed CMM-vs-GSM bound/stoichiometry swaps to enable biomass. |
| `reintroduced_for_biomass.tsv` | Blacklisted reactions that had to be re-introduced after the swap-search to enable biomass. |

## Workflow

1. Open `BuildFSP237_ATPSafe_GSM.ipynb` in Jupyter.
2. Cell 1: imports + path constants. Set `BASE` to your local
   `imm904CobraModel` checkout if running outside the original env.
3. Cells 2–9: load CMM, verify baseline ATP, identify GSM-unique reactions,
   run ATP-safe expansion (~100-reaction batches with bisection),
   add biomass reaction (`bio_gsm`), find minimum swap set, re-introduce
   from blacklist if needed, and save.

All heavy lifting is in `atp_safe_expand.py`:
- `atp_yield(model, atp_rxn, glc_uptake)` — ATP yield FBA.
- `expand_atp_safe(...)` — iterative ATP-safe addition with batch + bisection.
- `minimum_shared_swaps_for_biomass(...)` — greedy minimum swap finder for biomass.
- `minimum_blacklist_for_biomass(...)` — re-add blacklisted reactions if biomass = 0 even after swaps.

## What it produces

`fsp237_atp_safe_gsm.json` — a genome-scale FSP237 model with:
- Aerobic ATP/glc = 30 (matches CMM)
- Anaerobic ATP/glc = 2 (matches CMM)
- All GSM-unique reactions that don't break ATP yields
- The CMM `bio_gsm` biomass reaction copied over

This file is the **input** to:
- `../fsp237_biomass_extension/extend_biomass.py` (next stage)
- `../gpr-update/apply_full_gpr_mapping.py` (parallel post-processing)

## Dependencies

- Python 3.10+ with `cobra`, `cobrakbase`, `pandas`, `numpy`
- KBase token at `~/.kbase/token` (only needed for the initial CMM fetch
  — the cached `cmm_169876_170_5.json` lets you skip this if present)
