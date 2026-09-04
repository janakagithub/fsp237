# Phase 3 (kickoff) — kofam KO annotation → EC map

kofamscan run on the **full FSP237 proteome** (14,857 proteins) via the poplar
Celery queue (`kofam_scan.exec_annotation`, KO db `2025-11-03`, 40 threads).
Task `85152ca3` → **SUCCESS** (rc 0). Raw outputs on shared scratch:
`/scratch/bioseed_tools/kofam_scan/run/fsp237_fullproteome_kofam_20260903/`
(`output` = mapper, `tabular/tabular.txt` = full hmmsearch, 120 MB). Local copies
in `outputs/kofam_mapper.tsv`, `outputs/kofam_tabular.txt`.

## KO → EC

EC pulled deterministically from KEGG `list/ko` definitions (`[EC:...]`),
`outputs/kegg_ko_list.tsv` (28,453 KOs; 11,134 carry EC). No LLM.

## Result — `outputs/gene_ko_ec.tsv` (union of kofam + xlsx KO)

| metric | count |
|---|---|
| genes with KO (kofam ∪ xlsx) | 5,281 |
| genes with EC | 2,659 |
| kofam-only genes | 980 |
| xlsx-only genes | 238 |
| both agree/overlap | 4,063 |

kofam added **980 genes** with KO that the `annota_final3_new.xlsx` lacked.

## Model coverage (V10, 1,274 genes / 1,045 GPR reactions)

| metric | count | % |
|---|---|---|
| model genes with KO | 934 | 73.3% |
| model genes with EC | 744 | 58.4% |
| **GPR reactions with ≥1 EC-bearing gene** | **917 / 1,045** | **87.8%** |

The 87.8% is the kcat-parameterization ceiling from EC alone. Reactions with a
gene but no EC (and the ~128 uncovered GPR reactions) fall to the kcat cascade's
sequence-based tier (DLKcat/TurNuP) or EC-family/default fallback.

## Files written

- `outputs/kofam_mapper.tsv` — gene → KO (thresholded)
- `outputs/kofam_tabular.txt` — full hmmsearch table (scores/thresholds)
- `outputs/kegg_ko_list.tsv` — KEGG KO definitions (KO→EC source)
- `outputs/gene_ko_ec.tsv` — **canonical** gene → KO(s) → EC(s) + source
- `outputs/gene_annotation_map.tsv` — updated with `ko_all`, `ec_number`, `ec_ko_source`
- `outputs/kofam_job.json` — job provenance

## Next

kcat parameterization: EC → BRENDA/SABIO-RK kcat (EC hierarchy, organism-preferred
→ fungal → any), gap-fill with DLKcat/TurNuP on the proteome sequences, then
EC-family median / default fallback. Feeds the sMOMENT baseline, then GECKO.
