# `fsp237_biomass_extension/` — extend the yeast-derived biomass to fit FSP237

**Stage 2** of the model-build pipeline. The ATP-safe GSM from stage 1
inherits the iMM904 yeast biomass reaction, which is missing several
*Colletotrichum sublineola*–specific cell-wall and storage compounds.
This stage curates those in.

It also hosts the **site-build scripts** for the published static site at
https://janakagithub.github.io/fsp237/atp-safe/ and the **infection-stage
simulation plan** that drives the carbon-source panel.

## Files in this directory

### Biomass extension

| File | What it is |
|---|---|
| **`extend_biomass.py`** | The main script. Reads the curated Excel of missing reactions + a hardcoded list of fungal pathway additions (mannitol biosynthesis, α-1,3-glucan synthase, chitin route corrections), adds them to the model, then extends `bio_gsm` with the Colletotrichum-specific biomass demands (chitin, melanin, mannitol, α-1,3-glucan, glycogen rebalanced, dNTPs rebalanced for GC content). |
| `15June2026_Missing reactions GEM_higginsianum.xlsx` | Input — curated list of reactions to add from C. higginsianum literature. |
| `biomass_extension_log.tsv` | Output — every biomass-component change with reason. |
| `extension_report.tsv` | Output — every reaction added/blocked with status, equation, GPR. |
| **`fsp237_atp_safe_gsm_extended.json`** | Output model — ATP-safe GSM + curated additions + extended biomass. |
| `INFECTION_SIM_PLAN.md` / `.docx` | The simulation plan: 18 carbon-source conditions across pre-infection / biotrophic / necrotrophic / cocktail stages. Drives `../simulations/run_simulation_panel.py`. |

### Site builders

| File | What it builds |
|---|---|
| **`build_atp_safe_site.py`** | Generates `/atp-safe/reactions.json` (the data file the website fetches). Reads the current canonical V10 model, runs aerobic + anaerobic FBA, embeds reaction list, biomass composition, simulation panel results (joined from `../simulations/simulation_results.tsv`), and per-condition literature support (from `../simulations/condition_literature.tsv`). |
| **`build_escher_maps.py`** | Generates `/atp-safe/map_aerobic.html` and `/atp-safe/map_anaerobic.html` — the compartment-coloured Escher flux maps (cytosol pink, mito green, peroxisome orange, etc., intensity scaled by |flux|). |

Both scripts hardcode paths to the local checkout — adjust the `BASE`
constant at the top if running elsewhere.

## Workflow

```
extend_biomass.py
   ↓
fsp237_atp_safe_gsm_extended.json   →  ../gpr-update/apply_full_gpr_mapping.py
                                          ↓
                                       fsp237_atp_safe_gsm_gpr_updated.json
                                          ↓
                                       ../simulations/gapfill_v1_v2/build_v1..v10_*.py
                                          ↓
                                       fsp237_gapfilled_Version10_*.json
                                          ↓
                                       build_atp_safe_site.py + build_escher_maps.py
                                          ↓
                                       atp-safe/reactions.json + map_*.html → site
```

## What's in `bio_gsm` after this stage

74 components total. The Colletotrichum extensions (highlighted with the
`extension` badge on the site Biomass tab):

| Category | Compound | Coefficient | Citation |
|---|---|---:|---|
| Cell wall | Chitin (cpd11683_c0) | −0.05 | rxn15558 chitin synthase route |
| Cell wall | 1,3-α-D-Glucan (cpd12148_c0) | −0.05 | rxn15561 Ags1 (gene_5001); Fujikawa 2012 PMID 22927818 |
| Storage | Glycogen (cpd00155_c0) | −0.30 | rebalanced from KBase default 0.5185 |
| Storage | D-Mannitol (cpd00314_c0) | −0.20 | fungal MpdA pathway (rxn00546 + rxn01560) |
| Pigment | Melanin (cpd12744_c0) | −0.001 | DHN-melanin pathway (rxn00024 → rxn06809) |
| DNA | dAMP, dTMP | −0.00284 each | rebalanced to GC=52.7% per TX430BB (Baroncelli 2014) |
| DNA | dGMP, dCMP | −0.00316 each | rebalanced to GC=52.7% |

## Dependencies

- Python 3.10+ with `cobra`, `pandas`, `openpyxl`
- For the site builders: `escher` 1.8+ (only `build_escher_maps.py`)
