# Phase 1 — proteomics-driven E-Flux

Model fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json; E-Flux (Colijn 2009) with LINEAR TMT intensity, per-condition P99 normalization; GPR max_OR(min_AND). Baseline = unconstrained pFBA on the same medium.

## Per-condition summary

```
condition       media_key  glucose_uptake  baseline_biomass  eflux_biomass eflux_status  biomass_ratio  n_active_baseline  n_active_eflux  n_rxn_scaled  n_rxn_throttled   p99_linear
      PDA 19_pdb_baseline             5.0          0.191472       0.191472      optimal            1.0                509             509           882              869 1.361140e+10
     half     20_pdb_half             2.5          0.095736       0.095736      optimal            1.0                509             509           882              869 1.517621e+10
 onetenth 21_pdb_onetenth             0.5          0.019147       0.019147      optimal            1.0                507             507           882              869 8.648133e+09
```

## Capacity titration (proteome-implied bottleneck onset)

Biomass as a fraction of own-medium baseline as the global capacity multiplier C shrinks. First bottleneck reaction(s) per condition: {"PDA": ["rxn08617_c0"], "half": ["rxn00198_m0", "rxn08617_c0"], "onetenth": ["rxn00198_m0", "rxn08617_c0"]}

## E-Flux ↔ proteome DE concordance

Reactions with |Δflux|>0 and |proteome logFC|>0.5; fraction whose E-Flux flux change agrees in sign with the proteome fold-change.

- **half_vs_PDA**: 40/120 = 33.3% concordant
- **onetenth_vs_PDA**: 82/163 = 50.3% concordant
- **onetenth_vs_half**: 86/160 = 53.8% concordant

## Outputs
- `<cond>_eflux.tsv`, `eflux_flux_matrix.tsv`, `eflux_summary.tsv`, `eflux_de_concordance.tsv`