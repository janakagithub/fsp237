# FSP237 GPR overhaul — report

Generated 2026-06-18.

## Objective

Replace foreign gene identifiers (CH63R_* from C. higginsianum curation,
Y####W/C yeast IDs from iMM904, NP_* RefSeq accessions) in the FSP237 model's
GPRs with FSP237-native `gene_*` identifiers, using protein BLAST to find the
best ortholog.

User rules:
- `gene_*` tokens (FSP237 native) are preserved verbatim — never overwritten.
- CH63R_* tokens are replaced with the best `gene_*` BLAST hit.
- If a CH63R has no acceptable hit (pident ≥ 30, qcov ≥ 50, e-value ≤ 1e-10),
  it is removed and the reaction is flagged for manual review.
- MSTRG (novel-transcript) hits are reported in the table but NOT added to
  the model — only `gene_*` proteins (12,527 of 14,857 in the proteome FASTA)
  are eligible as final mapping targets.
- All other foreign tokens (yeast, NP_) are stripped.

## Pipeline

1. **Pull the FSP237 KBase genome** (`172651/44/3`) — gives `gene_NNNN` IDs
   with their protein translations (the missing link between the model and
   the gbff).
2. **Extract CH63R proteins** from `C_higgensium.gbff` (564 unique CH63R IDs
   across both curated Excels; 563 had sequences).
3. **Build BLAST DB** from
   `Csublineola_reference_plus_novel_classu.proteins.fa`
   (12,527 gene_* + 2,330 MSTRG-novel = 14,857 sequences).
4. **Run blastp** for each CH63R query against the full DB; capture identity,
   coverage, e-value, bitscore.
5. **Pick best gene_* per CH63R** subject to thresholds.
6. **Map curated CH63R sets per reaction** from both Excel files:
   - `Integration_Higgginsianum_GEM - 4JK.xlsx`
   - `15June2026_Missing reactions GEM_higginsianum.xlsx`
7. **Update model GPRs**: union of preserved `gene_*` + mapped CH63R, deduped.
8. **Push updated model** to the visualizer.

## BLAST results

| | count |
|---|---|
| CH63R queries with sequence | 556 |
| Mapped to `gene_*` above thresholds | **547 (98.4%)** |
| No hit passing thresholds (flagged) | 9 |
| Weak matches kept (pident < 50%, flagged) | 34 |
| MSTRG was best-overall (gene_* still used) | 56 |

## Model GPR changes

| | count |
|---|---|
| Reactions updated | **514** |
| Reactions flagged (issues during mapping) | 99 |
| Total CH63R → `gene_*` substitutions | 981 |
| Foreign yeast/NP tokens removed | 143 |
| Unmapped CH63R dropped from GPRs | 34 |

After the update:
- `0` reactions contain CH63R tokens
- `0` reactions contain yeast Y####W/C or NP_* tokens
- Model gene count: 1128 → 1286 (158 newly-referenced FSP237 genes from BLAST)

## Validation (preserved)

| | aerobic | anaerobic | target |
|---|---|---|---|
| ATP yield (glc=1) | 30.000 | 2.000 | 30 / 2 |
| Biomass (GSM-style media) | 0.1809 | 0.0369 | ~0.181 / ~0.037 |

## Deliverables (in `/gpr-update/`)

- `fsp237_atp_safe_gsm_gpr_updated.json` — updated model (use this going forward)
- `ch63r_to_fsp237_mapping_full.tsv` — full mapping table with BLAST evidence
- `gpr_change_log.tsv` — per-reaction before/after GPR change record
- `flagged_reactions.tsv` — 99 reactions needing manual review
- `blast_results_all.tsv` — raw blastp hits (1,871 lines)
- `ch63r_queries.faa` / `ch63r_queries_v2.faa` — query FASTAs
- `blast_db/` — BLAST protein DB of FSP237 + novel
- `apply_full_gpr_mapping.py` — driving script (idempotent: re-run anytime)
- `apply_gpr_mapping.py` — initial smaller script (kept for reference)

## Visualizer

Live at https://janakagithub.github.io/fsp237/atp-safe/

GPRs in the Reactions and Active fluxes tabs now show:
- FSP237 `gene_NNNN` identifiers throughout
- Long "or" lists wrap cleanly within the GPR column (max-width 260 px)
- EC numbers (newly added by reading from reaction annotations when not in
  the static TSV)

## Next step

When ready, save the updated model to KBase. Per user instruction, hold off
until the biomass extension work (chitin/melanin/mannitol/etc.) is finalized.
