# `gpr-update/` — replace foreign gene identifiers with FSP237 native `gene_*` IDs

**Stage 3** (parallel post-processing of `fsp237_biomass_extension/`).
The biomass-extended model still carries foreign gene identifiers in its
GPRs from upstream curation:

- `CH63R_*` from the *C. higginsianum* Excel curation
- `Y####W/C` from the iMM904 yeast template
- `NP_*` RefSeq accessions
- A few others

This pipeline replaces all of them with FSP237 native `gene_NNNN`
identifiers via protein BLAST against the FSP237 proteome.

## Files in this directory

### Scripts

| File | What it is |
|---|---|
| **`apply_full_gpr_mapping.py`** | The driver. Builds a BLAST CH63R→gene_* mapping from `blast_results_all.tsv`, reads the per-reaction CH63R gene sets from both curated Excels, walks every model reaction, replaces foreign GPR tokens with the best `gene_*` ortholog above quality thresholds (pident≥30, qcov≥50, e-value≤1e-10), and saves the updated model + change logs. |
| `apply_gpr_mapping.py` | Earlier prototype — kept for reference, not used in the canonical workflow. |

### BLAST inputs/outputs

| File | What it is |
|---|---|
| `ch63r_queries.faa` | First-pass FASTA of CH63R query proteins (initial set extracted from the iMM904-derived model). |
| `ch63r_queries_v2.faa` | Final FASTA of all 556 CH63R queries with sequences, pulled from `C_higgensium.gbff`. |
| `blast_db/` | BLAST protein DB of the FSP237 proteome (12,527 `gene_*` + 2,330 MSTRG-novel = 14,857 sequences). Pre-built; ship as-is so BLAST runs without prep. |
| `blast_results.tsv` | First-pass BLAST output (initial query set). |
| `blast_results_all.tsv` | Final BLAST output (1,871 lines covering all 556 queries against the FSP237 DB). |
| `blast_results_v2.tsv` | Intermediate refinement (kept for traceability). |
| `mstrg_hits.tsv` | CH63Rs whose best overall hit was an MSTRG (novel-transcript) protein. Not used as final mappings but recorded for traceability. |

### Mapping + output

| File | What it is |
|---|---|
| `ch63r_to_fsp237_mapping.tsv` | First-pass per-CH63R best-gene mapping table. |
| **`ch63r_to_fsp237_mapping_full.tsv`** | Final per-CH63R best-gene mapping table with all BLAST evidence (pident, qcov, e-value, bitscore) + flag column. |
| **`fsp237_atp_safe_gsm_gpr_updated.json`** | The updated model. 1703 reactions, 1286 genes (was 1128 — 158 new gene_* added via BLAST), zero foreign tokens. |
| `gpr_change_log.tsv` | Per-reaction before/after GPR change record (514 reactions updated, 981 CH63R substitutions, 143 yeast/NP tokens removed). |
| `flagged_reactions.tsv` | 99 reactions with issues during mapping (weak BLAST hits, unmapped CH63R dropped, GPR became empty). |
| `REPORT.md` | Narrative summary of the BLAST mapping outcome. |

## What's NOT in this directory (must be downloaded separately)

Two genome-annotation files were too large to commit (>50 MB combined
over GitHub's recommended single-file size). To re-run the BLAST yourself,
download them:

| File | Size | Source |
|---|---|---|
| `C_higgensium.gbff` | 85 MB | NCBI RefSeq `GCF_001672515.1` (C. higginsianum IMI 349063 RefSeq annotation) |
| `KBase_derived_Colletotrichum_sublineola_FSP237_MyCoCosm.gbff` | 116 MB | KBase / JGI MyCoCosm export |
| `Csublineola_reference_plus_novel_classu.proteins.fa` | 8 MB | (FSP237 proteome — used to build the BLAST DB; the DB itself IS shipped in `blast_db/`, so this isn't needed for re-running the mapping, only for re-building the DB) |

Once `C_higgensium.gbff` is in place, the script will work end-to-end. If
you only want to apply the existing mapping (use `blast_results_all.tsv` +
`fsp237_atp_safe_gsm_extended.json` from the previous stage), the script
runs without re-doing BLAST.

## Workflow

```
fsp237_atp_safe_gsm_extended.json  +  C_higgensium.gbff (download)
        │
        ▼  apply_full_gpr_mapping.py
        │  (Step 1: build CH63R -> gene_* mapping from BLAST)
        │  (Step 2: per-rxn CH63R index from both Excels)
        │  (Step 3: update each rxn's GPR — keep gene_*, map CH63R, drop foreign)
        ▼
fsp237_atp_safe_gsm_gpr_updated.json   ← canonical for downstream stages
ch63r_to_fsp237_mapping_full.tsv
gpr_change_log.tsv
flagged_reactions.tsv
```

## Numbers at the end of this stage

- 1703 reactions (unchanged)
- 1286 genes (was 1128 — 158 newly mapped from BLAST)
- 0 foreign tokens in any GPR
- 514 reactions updated; 981 CH63R substitutions; 143 yeast/NP tokens removed; 34 unmappable CH63R dropped
- Biomass + ATP yield invariants preserved (0.181 aer / 0.037 ana; 30 / 2)

## What happens next

`fsp237_atp_safe_gsm_gpr_updated.json` is the input to the gap-fill
pipeline in `../simulations/gapfill_v1_v2/`.
