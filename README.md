# FSP237 — *Colletotrichum sublineola* genome-scale metabolic model

A genome-scale metabolic model (GEM) of *Colletotrichum sublineola* strain
FSP237, the causal agent of sorghum anthracnose, with an infection-stage
simulation panel and a versioned gap-fill chain.

- **Live site (browse model + simulations interactively):** https://janakagithub.github.io/fsp237/atp-safe/
- **KBase model:** `28277/287/1` — `fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated` (linked to *C. sublineola* JGI genome `169876/166/3`)
- **Genome proxy for parameters:** *C. sublineola* TX430BB (Baroncelli et al. 2014 *Genome Announc* PMID 24926053)

## Repository layout

| Directory | Stage | What's there |
|---|---|---|
| **[`fsp237_atp_safe_gsm/`](fsp237_atp_safe_gsm/)** | 1 | Build the ATP-safe GSM from a curated central-carbon model (CMM) without breaking the 30/2 ATP yields. Produces `fsp237_atp_safe_gsm.json`. |
| **[`fsp237_biomass_extension/`](fsp237_biomass_extension/)** | 2 | Extend the yeast-derived biomass with Colletotrichum-specific compounds (chitin, melanin, mannitol, α-1,3-glucan, dNTP GC-rebalance, glycogen update). Hosts the **site builders** (`build_atp_safe_site.py` + `build_escher_maps.py`) and the **infection-stage simulation plan** (`INFECTION_SIM_PLAN.md`). |
| **[`fsp237_minimal_glucose/`](fsp237_minimal_glucose/)** | 0 (predecessor) | Initial glucose-viable build from raw KBase model, plus the ETC + futile-cycle repair scripts. Foundational notebook + companion fix-scripts. |
| **[`gpr-update/`](gpr-update/)** | 3 | Replace foreign gene IDs (CH63R / yeast Y####W/C / NP_*) with FSP237 native `gene_*` IDs via protein BLAST. Ships the pre-built BLAST DB. Produces `fsp237_atp_safe_gsm_gpr_updated.json`. |
| **[`simulations/`](simulations/)** | 4 | 18-condition × 2-O₂ FBA panel across pre-infection / biotrophic / necrotrophic / cocktail stages. Drives the gap-fill chain V1 → V10. **Current canonical model lives at `simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version10_*.json`.** |
| `atp-safe/` | — | The rendered static site (auto-built from V10 by `build_atp_safe_site.py` + `build_escher_maps.py`). |
| `CHECKPOINT.md` (in the source dir, not the repo) | — | Snapshot of every reference (GitHub SHA, KBase ref, file paths) for resuming work — see local checkout. |

## End-to-end build chain

```
KBase 169876/251/1 ── BuildMinimalFSP237Model.ipynb ── fsp237_minimal_glucose.json
  (raw FSP237)            (+ ETC + futile-cycle fixes)        │
                                                              │ (media bounds reference)
                                                              ▼
KBase 169876/170/5 ── BuildFSP237_ATPSafe_GSM.ipynb ── fsp237_atp_safe_gsm.json
  (curated CMM)           (ATP-safe expansion)                │
                                                              ▼
                          extend_biomass.py ── fsp237_atp_safe_gsm_extended.json
                                                              │
                                                              ▼
                  apply_full_gpr_mapping.py ── fsp237_atp_safe_gsm_gpr_updated.json
                                                              │
                                                              ▼
  build_v1_gapfill.py → V1 → build_v2_integrate_genes.py → V2
                                                              │
                              build_v3v4_dedup.py → V3 / V4   │
                                                              │
                              build_v5v6_dirlock.py → V5 / V6
                                                              │
                              build_v9v10_vlcfa_cleanup.py → V9 / V10  ← current
                                                              │
                                                              ▼
  build_atp_safe_site.py + build_escher_maps.py → atp-safe/{reactions.json, map_*.html}
                                                              │
                                                              ▼
                                                        Live website
                                                              │
                                                              ▼
                                                       save_v6_to_kbase.py → KBase 28277/287/1
```

Each stage produces a JSON model that the next stage consumes. Every
intermediate is preserved on disk (and most are on GitHub) so any step
can be re-run or audited independently.

## Quick start

```bash
git clone https://github.com/janakagithub/fsp237.git
cd fsp237
```

To load the **current canonical model** in COBRApy:

```python
import cobra
m = cobra.io.load_json_model(
    'simulations/gapfill_v1_v2/models/'
    'fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json'
)
print(len(m.reactions), 'reactions /', len(m.metabolites), 'metabolites /', len(m.genes), 'genes')
# 1622 reactions / 1268 metabolites / 1274 genes
```

To **re-run the 18-condition simulation panel** against V10:

```bash
python simulations/run_simulation_panel.py
# writes simulations/simulation_results.tsv + per-condition flux dumps
```

To **read the master narrative** of every gap-fill decision:

```
simulations/gapfill_v1_v2/reports/SUMMARY.md
```

To **see the pathway diagrams** (β-ox, Penttilä L-Ara, Ashwell D-galU,
oleate):

```
simulations/gapfill_v1_v2/reports/PATHWAY_DIAGRAMS.md      (+ .docx)
simulations/gapfill_v1_v2/reports/OLEATE_PATHWAY.md         (+ .docx)
```

## Dependencies

- Python 3.10+
- `cobra`, `pandas`, `numpy`, `openpyxl`
- For site builders: `escher` 1.8+
- For BLAST: `blastp` on PATH (NCBI BLAST+)
- For KBase save: `cobrakbase`, KBase token at `~/.kbase/token`

## Large files not in this repo

Two genome-annotation files are too large for GitHub. To re-run the full
GPR-update pipeline (rather than just applying the existing mapping):

| File | Size | Source |
|---|---|---|
| `gpr-update/C_higgensium.gbff` | 85 MB | NCBI RefSeq `GCF_001672515.1` |
| `gpr-update/KBase_derived_Colletotrichum_sublineola_FSP237_MyCoCosm.gbff` | 116 MB | KBase / JGI MyCoCosm |

The BLAST DB (`gpr-update/blast_db/`) **is** shipped, so you can re-do
mapping decisions or BLAST new queries without rebuilding the DB.

## Citation

The major model build, GPR overhaul, biomass extension, and infection-stage
simulation panel are documented in:

- `simulations/gapfill_v1_v2/reports/SUMMARY.md` — full version chain
- `fsp237_biomass_extension/INFECTION_SIM_PLAN.md` — simulation plan + biology
- `simulations/condition_literature.tsv` — 40 PMID citations per condition

Underlying genome reference: Baroncelli R, Sanz-Martín JM, Rech GE, Sukno SA,
Thon MR. 2014. *Draft Genome Sequence of Colletotrichum sublineola, a
Destructive Pathogen of Cultivated Sorghum.* Genome Announc 2(3):e00540-14.
**PMID 24926053**.
