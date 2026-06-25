# FSP237 ATP-safe GSM expansion — final report

Generated 2026-06-16.

## Strategy

Start with the KBase-curated central carbon model (CMM) `169876/170/5` for *C. sublineola* FSP237,
which produces **exactly 30 mmol ATP / mmol glucose aerobically and 2 anaerobically**. Iteratively
add reactions from the larger fsp237 GSM (`fsp237_minimal_glucose.json`, post-ETC fixes), rejecting
(blacklisting) any reaction that pushes ATP yield above the curated values — the signature of an
introduced futile cycle.

**CMM has precedence**: shared reactions (same ID in both models but different bounds/stoichiometry)
keep the CMM's version by default. They only get swapped to the GSM's version if biomass cannot
grow otherwise.

## Stage 1: ATP-safe iterative expansion

CMM and GSM share 44 reaction IDs. Of the 1670 GSM reactions, 1474 are GSM-unique (excluding `bio1`).
Each candidate was tested by batched addition (size 100) and bisected on any batch that pushed the
(aerobic, anaerobic) ATP yield above (30+1e-3, 2+1e-3).

| | count |
|---|---|
| Candidate reactions (GSM-unique, excl. bio1) | 1474 |
| Safely added (ATP yield preserved) | **1468** |
| Blacklisted (would have exceeded ATP yield) | **6** |

### Blacklisted reactions (6)

| rxn_id | name | why |
|---|---|---|
| `rxn08096_c0` | α-ketoglutarate diffusion (c0↔x0) | reversible transporter creates a peroxisomal aKG loop |
| `rxn00903_m0` | L-Val:2-oxoglutarate aminotransferase | reversible without coupled NADH; enables a transamination cycle |
| `EX_cpd00159_e0` | L-lactate exchange | secretion enables a lactate-fermentation alt path that perturbs anaerobic yield |
| `rxn08016_c0` | palmitate-ACP ligase | bacterial-style acyl-ACP synth, no fungal homolog; reversible bounds form ATP cycle |
| `rxn34741_c0` | long-chain FA-ACP ligase (14:0) | same pattern as rxn08016_c0 |
| `rxn34744_c0` | long-chain FA-ACP ligase (18:1) | same pattern as rxn08016_c0 |

All six are bacterial-style acyl-ACP / acyl-CoA ligases or unannotated transporters / transaminases
running reversibly — same dedup-miss pattern as the previously-fixed `rxn08017_c0` / `frxn08975_c0`
duplicates. See `blacklist.tsv` for full equations.

## Stage 2: biomass test

The GSM biomass reaction (70 metabolites) was added as `bio_gsm` in the expanded model and tested
on GSM-style minimal media (glucose -5, O2 -10, sulfate/Fe/Na/K/NH3/Pi/H2O/H⁺/CO2 open).

**Biomass without any shared-rxn swaps: 0** — at least one CMM reaction is blocking growth.

## Stage 3: minimum shared-rxn swap for biomass

Algorithm: try all 30 shared-rxn-mismatches swapped to GSM versions; if biomass grows, greedily
revert each swap and keep only the ones whose reversion drops biomass to zero.

**Result: only 1 shared reaction needed to swap.**

| rxn_id | issue | resolution |
|---|---|---|
| `HACNHm_m0` | CMM uses metabolite `Mhicitm_m0` (a CMM-local id for homoisocitrate); GSM uses `hicit-m_m0` (the canonical id used everywhere else in the GSM, including by `HICITDm_m0` and the lysine biosynth pathway). The CMM's `Mhicitm_m0` is an orphan — no other reaction in the merged model produces or consumes it. | Use GSM version so `HACNHm_m0` connects to the rest of the lysine pathway. |

The other 29 shared reactions with bounds/stoich mismatches kept their CMM versions. Notably,
`rxn35242_c0` (pyruvate kinase) and `rxn30381_c0` (glycerol-3-P dehydrogenase) — which initially
*looked* like they should be swapped — actually work as-is because the CMM has them written with
reversed metabolite-side ordering AND reverse-only bounds, so they run in the biologically-correct
direction (PEP+ADP→Pyr+ATP at flux −6.56 in the final solution).

## Stage 4: final verification

| condition | value | target |
|---|---|---|
| Aerobic ATP/glucose | **30.000** | 30 |
| Anaerobic ATP/glucose | **2.000** | 2 |
| Aerobic biomass flux | **0.1829** | > 0 |
| Anaerobic biomass flux | **0.0371** | > 0 |

The expanded model preserves the CMM's curated ATP yield exactly **AND** grows biomass aerobically
and anaerobically on standard minimal-glucose media, **without re-introducing any blacklisted
reaction**.

The aerobic biomass (0.183) is ~19% lower than the unconstrained GSM's biomass (0.226). That gap
reflects what was previously "free ATP" from the 6 blacklisted futile cycles silently inflating
growth in the original GSM.

## Reaction accounting

| | count |
|---|---|
| Reactions in expanded model | **1688** |
| Metabolites in expanded model | **1317** |
| Safely added from GSM | 1468 |
| Blacklisted (excluded) | 6 |
| Re-introduced from blacklist for biomass | 0 |
| Shared rxns swapped CMM→GSM (minimum needed) | 1 |
| Shared rxns kept as CMM (precedence respected) | 29 |
| CMM-unique reactions (kept as-is) | 175 |

## Files

- `fsp237_atp_safe_gsm.json` — the expanded, ATP-safe GSM (use this for downstream analysis)
- `cmm_169876_170_5.json` — cached CMM as fetched from KBase 169876/170/5
- `safe_added.tsv` — 1468 safely-added reactions
- `blacklist.tsv` — 6 blacklisted reactions
- `shared_rxns_swapped.tsv` — 1 shared reaction where the GSM version was needed
- `reintroduced_for_biomass.tsv` — empty (no re-introduction needed)
- `BuildFSP237_ATPSafe_GSM.ipynb` — driving notebook
- `atp_safe_expand.py` — helper module

## Conclusion

Successfully built an expanded GSM (1688 reactions, 1317 metabolites) that:
- Preserves the CMM's energy-accurate ATP yield (30 aerobic / 2 anaerobic per glucose), exactly.
- Grows biomass aerobically (0.183) and anaerobically (0.037) on standard minimal-glucose media.
- Respects CMM precedence — only 1 of 30 shared-reaction mismatches required swapping to GSM.
- Has 6 ATP-cycle reactions documented in the blacklist (not silently inflating ATP yield).
- Required 0 blacklist re-introduction — the model is fully biomass-feasible while preserving the
  curated energy phenotype.
