# FSP237 project — checkpoint for resuming work later

**Snapshot taken (latest):** 2026-07-26
**Project**: FSP237 (*Colletotrichum sublineola*) genome-scale metabolic model
for studying anthracnose infection on sorghum.

---

## The single most important identifier

**GitHub commit SHA: `007b0a7…`** (short form: **`007b0a7`**, 2026-07-26)
— "Stage 5 pathway analysis + expression-painted Escher map + Pathway
breakdown tab". Preceded by NaN-fix `5071fc9` and Stages 0/1/2/4 landing
commit `1ae0085` (2026-07-24).

This commit on https://github.com/janakagithub/fsp237 (branch `main`) is
the **canonical reference for "where the FSP237 work is right now"**. If
you ever need to pin or revert, that's the SHA.

To restore the full repo state at this checkpoint:
```bash
git clone https://github.com/janakagithub/fsp237.git
cd fsp237
git checkout 1ae00851d4ed48e07afb6b88355c6f746b18987d
```

**Commit chain (most recent first)** — `git log --oneline`:
- `007b0a7` — Stage 5 pathway analysis + expression-painted Escher map + Pathway breakdown tab (2026-07-26)
- `5071fc9` — fix(site): strip NaN from reactions.json — was breaking JSON.parse in browsers
- `3a99c1b` — atp-safe: RNA-seq tab polish — top-fit callout + GitHub raw-output links
- `058be85` — CHECKPOINT.md: bump to 1ae0085
- `1ae0085` — rna_seq_integration: S1 transcriptomics × V10 GEM (Stages 0/1/2/4) + new RNA-seq tab (2026-07-24)
- `a23ebfe` — simulations: update condition_literature.tsv to Cs-first set (45 cites + relevance_tier column)
- `51606df` — atp-safe: Cs-first literature audit + relevance tier badges
- `e40cd6b` — Add dedup_initial_build.py (extracted from BuildMinimalFSP237Model.ipynb)
- `596a64b` — Add full FSP237 model-build pipeline + simulations (389 files, ~36 MB; all major scripts uploaded with per-folder READMEs)
- `d678534` — atp-safe: surface biomass-extension compounds on Biomass tab
- `89e56e1` — atp-safe: add Literature column to Simulations tab (40 PubMed-linked citations)
- `06090c2` — atp-safe: pivot Simulations tab -- aerobic + anaerobic side-by-side
- `0a6961b` — atp-safe: publish V10 -- complete VLCFA chain + unused-compartment cleanup (previous checkpoint)

The repo now contains both the rendered static site (`atp-safe/`) AND the
full source pipeline (`fsp237_atp_safe_gsm/`, `fsp237_biomass_extension/`,
`fsp237_minimal_glucose/`, `gpr-update/`, `simulations/`) with READMEs.

---

## Where the work lives (4 stable surfaces)

### 1. GitHub repo — the published static site + commit history

- **Repo**: https://github.com/janakagithub/fsp237 (main branch)
- **Live site**: https://janakagithub.github.io/fsp237/atp-safe/
- **Latest commit at checkpoint**: `0a6961b` (V10 publish)
- **Files tracked**: `atp-safe/index.html`, `atp-safe/reactions.json`,
  `atp-safe/map_aerobic.html`, `atp-safe/map_anaerobic.html`, plus older
  root-level `index.html`, `map.html`, `reactions.json`
- **Workflow**: every site refresh is a commit with a descriptive message
  (mannitol pathway, GPR overhaul, alpha-1,3-glucan, V10 VLCFA, etc.) — so
  the commit log is the project history. `git log --oneline` reads as a
  changelog.

### 2. KBase workspace 28277 — model snapshots

Workspace: `janakakbase:narrative_1518190880851` ("Fungal Biomass testing
Narrative"). Two FSP237 models from this checkpoint:

| Object ref | Name | What it is |
|---|---|---|
| **`28277/287/1`** | `fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated` | **Current — V10**, the one the site publishes. 1622 rxns / 1268 mets / 1274 genes. Full VLCFA chain, slim compartments, all GPRs intact. |
| `28277/286/3` | `fsp237_gapfilled_Version6_dirlock_genes_integrated` | Earlier V6 snapshot — pre-VLCFA, with g0/n0/v0 still in. Kept for comparison. |

Genome reference for both: **`169876/166/3`** (*C. sublineola* JGI, 13047
features — provides the gene_NNNN feature refs so KBase viewer renders
GPRs correctly).

Viewer URLs:
- V10: https://narrative.kbase.us/#dataview/28277/287/1
- V6 : https://narrative.kbase.us/#dataview/28277/286/3

### 3. Local working directory

**Root**: `/home/janakae/fungalTemplate/imm904CobraModel/`

Layout (key subdirectories):
```
imm904CobraModel/
├── CHECKPOINT.md                  ← this file
├── fsp237_minimal_glucose/        ← minimal-glucose GSM build + first flux map
│   └── BuildMinimalFSP237Model.ipynb
├── fsp237_atp_safe_gsm/           ← ATP-safe CMM-anchored GSM build
│   └── BuildFSP237_ATPSafe_GSM.ipynb
├── fsp237_biomass_extension/      ← chitin/melanin/mannitol/α-glucan biomass + Excel pipeline
│   ├── extend_biomass.py
│   ├── build_atp_safe_site.py     ← site reactions.json builder
│   ├── build_escher_maps.py       ← Escher map builder
│   ├── INFECTION_SIM_PLAN.md      ← Cs infection-stage simulation plan
│   └── INFECTION_SIM_PLAN.docx
├── gpr-update/                    ← CH63R → FSP237 BLAST mapping pipeline
│   ├── REPORT.md
│   ├── apply_full_gpr_mapping.py
│   ├── blast_db/                  ← FSP237 protein BLAST DB (reusable)
│   ├── Csublineola_reference_plus_novel_classu.proteins.fa
│   └── C_higgensium.gbff          ← CH63R gene annotations + sequences
├── simulations/
│   ├── run_simulation_panel.py    ← 18-condition × 2-O2 simulation driver
│   ├── simulation_results.tsv     ← latest panel results (V10)
│   ├── per_condition/             ← per-condition non-zero flux dumps
│   ├── RESULTS.md                 ← panel interpretation + gap-fill recipe
│   └── gapfill_v1_v2/             ← gap-fill V1→V10 + reports
│       ├── build_v1_gapfill.py    ← V1: 28 gap-fill rxns (no GPRs)
│       ├── build_v2_integrate_genes.py  ← V2: BLAST gene picks
│       ├── build_v3v4_dedup.py    ← V3/V4: exact-duplicate cleanup
│       ├── build_v5v6_dirlock.py  ← V5/V6: degradation-only direction lock
│       ├── build_v9v10_vlcfa_cleanup.py ← V9/V10: VLCFA + compartment prune
│       ├── save_v6_to_kbase.py    ← KBase push (parameterized for any V*)
│       ├── test_gapfilled_model.py← re-run panel against any model JSON
│       ├── find_candidates.py     ← C. higginsianum candidate finder + BLAST
│       ├── models/                ← V1, V2, V3, V4, V5, V6, V7, V8, V9, V10 JSON
│       ├── candidates/            ← BLAST query FASTAs + provenance
│       ├── blast/                 ← raw BLAST results + rxn→gene mapping
│       └── reports/
│           ├── SUMMARY.md         ← master narrative of the V1→V10 chain
│           ├── PATHWAY_DIAGRAMS.md / .docx  ← β-ox + L-Ara + Ashwell diagrams
│           ├── OLEATE_PATHWAY.md / .docx
│           ├── OLEATE_LITERATURE_SUPPORT.md / .docx
│           ├── v{1..10}_simulation_results.tsv  ← panel per version
│           ├── v{3,4,5}_dedup_log.tsv / _direction_locks.tsv
│           ├── v9_change_log.tsv
│           └── v{1..10}_per_condition/  ← per-condition flux dumps per version
└── (~445 files total in the working dir, including older intermediate work)
```

**To zip and archive this entire local state**:
```bash
cd /home/janakae/fungalTemplate
tar czf fsp237_checkpoint_2026-06-24.tar.gz imm904CobraModel/
# also include the publish repo if you want the rendered site:
tar czf fsp237_publish_repo_2026-06-24.tar.gz /home/janakae/fsp237/
```

### 4. Claude auto-memory (Claude Code sessions only)

Path: `/home/janakae/.claude/projects/-home-janakae-fungalTemplate-imm904CobraModel/memory/`

This is Claude Code's persistent memory across sessions in *this working
directory*. The files there are plain markdown and survive across Claude
sessions automatically. Not necessary for reproducing the work (everything
biological is in #1–#3), but useful so any future Claude session in this
directory picks up the context.

Files at checkpoint:
- `MEMORY.md` — index
- `project_fsp237_naming.md` — naming conventions
- `project_fsp237_vs_tx430bb.md` — genome proxy citation
- `project_gmm_glucose_cap.md` — glucose uptake bounds
- `project_fsp237_etc_state.md` — ETC fixes documented
- `reference_fsp237_publish.md` — full publish workflow + refresh order

These will continue to load when you return to this directory in a future
Claude Code session. They will NOT auto-load if you `cd` to a different
project. If you want to seed a different Claude project with FSP237
context, copy the relevant `.md` files into that project's
`.claude/projects/.../memory/` folder.

---

## Model version trail (full lineage, all preserved on disk)

| Version | File (under `simulations/gapfill_v1_v2/models/`) | What it adds |
|---|---|---|
| source | `gpr-update/fsp237_atp_safe_gsm_gpr_updated.json` | GPR-clean biomass-extended baseline |
| V1 | `fsp237_gapfilled_Version1_noGenes.json` | 28 gap-fill rxns (β-ox C8/C6/C4, cofactor shuttles, Penttilä L-Ara, Ashwell galU); no GPRs |
| V2 | `fsp237_gapfilled_Version2_gapfill_genes_integrated.json` | V1 + BLAST-derived GPRs on 19 rxns |
| V3 | `fsp237_gapfilled_Version3_dedup_noGenes.json` | V1 + exact-duplicate cleanup (34 rxns dropped) |
| V4 | `fsp237_gapfilled_Version4_dedup_genes_integrated.json` | V2 + same dedup, GPRs merged |
| V5 | `fsp237_gapfilled_Version5_dirlock_noGenes.json` | V3 + β-ox + Ashwell direction-locked to degradation only |
| V6 | `fsp237_gapfilled_Version6_dirlock_genes_integrated.json` | V4 + same dirlock — first KBase save (`28277/286/3`) |
| V7 | `fsp237_gapfilled_Version7_vlcfa_noGenes.json` | V5 + first (buggy) VLCFA attempt — kept for record |
| V8 | `fsp237_gapfilled_Version8_vlcfa_genes_integrated.json` | V6 + same (buggy) VLCFA — kept for record |
| **V10** | `fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json` | **Current** — V6 + complete VLCFA chain (C26→C24→C22→C20→C18, biosynthesis-consistent cpd IDs) + 88 unused-compartment rxns dropped (g0, n0, v0) — published on the site, saved to KBase as `28277/287/1` |
| V9 | `fsp237_gapfilled_Version9_vlcfa_complete_noGenes.json` | V5 + same V10 changes minus the gene assignments |

Each version is a self-contained COBRApy JSON model — load with
`cobra.io.load_json_model(path)`.

---

## How to resume work later

1. **Pull the repo + checkout the checkpoint commit** so the published-
   site state is reproducible:
   ```bash
   git clone https://github.com/janakagithub/fsp237.git
   cd fsp237 && git log -1   # confirm you're at 0a6961b or rebase fwd
   ```
2. **Restore the local working directory** from your archive (or `cp -r`
   from this checkpoint if still on disk):
   ```bash
   tar xzf fsp237_checkpoint_2026-06-24.tar.gz
   cd imm904CobraModel
   ```
3. **Open the master narrative** at
   `simulations/gapfill_v1_v2/reports/SUMMARY.md` — it walks through the
   V1→V10 chain with every decision, plus links to the per-version
   flux/dedup/dirlock/VLCFA logs.
4. **Load the current model**:
   ```python
   import cobra
   m = cobra.io.load_json_model(
       'simulations/gapfill_v1_v2/models/'
       'fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json'
   )
   ```
5. **Re-run the simulation panel** to verify the local state matches
   the published / KBase state:
   ```bash
   python simulations/run_simulation_panel.py
   # diff simulations/simulation_results.tsv against the version that's
   # embedded in /home/janakae/fsp237/atp-safe/reactions.json (under the
   # 'simulations' key)
   ```
6. **For Claude Code sessions** in a different directory, copy
   `~/.claude/projects/-home-janakae-fungalTemplate-imm904CobraModel/memory/*.md`
   into the new project's memory folder if you want the context to follow.

---

## Quick references (you can hand any of these to a colleague or to a future you)

| Use case | Where to point |
|---|---|
| Live site to browse | https://janakagithub.github.io/fsp237/atp-safe/ |
| Source-of-truth git commit | `0a6961bf5626d09459a6b3f0f8116860998af3cc` |
| KBase model for collaborators | `28277/287/1` (V10) |
| KBase narrative | https://narrative.kbase.us/narrative/1518190880851 |
| Current model JSON | `simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json` |
| Full project history | `git log` on `github.com/janakagithub/fsp237` + `SUMMARY.md` in this repo |
| Simulation panel + biological rationale | `simulations/gapfill_v1_v2/reports/SUMMARY.md` + `fsp237_biomass_extension/INFECTION_SIM_PLAN.md` |
| Pathway-by-pathway diagrams | `simulations/gapfill_v1_v2/reports/PATHWAY_DIAGRAMS.md` (also `.docx`) |
| Oleate-specific lit support | `simulations/gapfill_v1_v2/reports/OLEATE_LITERATURE_SUPPORT.md` |

---

## TL;DR

If you remember only one thing: **`github.com/janakagithub/fsp237` commit
`0a6961b`** and **KBase `28277/287/1`**. Those two refs let you (or
anyone else) reproduce the model, the website, and every simulation
result. Everything else in this repo is supporting material that can be
regenerated from those two anchors.
