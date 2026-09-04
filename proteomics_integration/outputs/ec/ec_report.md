# Phase 2 — proteome-allocated enzyme-capacity model

Model fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json. Catalyzed reactions costed: **1045**. **kcat from BRENDA 2025_1** (experimental turnover numbers, organism-tiered Colletotrichum→fungal→any; EC-family and Bar-Even EC-class prior as fallback). BRENDA-backed: **913/1045** (87.4%). Source breakdown — brenda_any: 406, brenda_fungal: 369, brenda_ecfamily: 137, default_no_ec: 128, ecclass_prior: 4, brenda_colletotrichum: 1. MW from FSP237 proteome sequences (median 44.1 kDa).

## Enzyme-budget growth curve (frac of own-medium baseline)

```
condition        PDA    half  onetenth
budget                                
1.000000e+10  1.0000  1.0000    1.0000
3.000000e+09  1.0000  1.0000    1.0000
1.000000e+09  1.0000  1.0000    1.0000
3.000000e+08  0.3973  0.6619    1.0000
1.000000e+08  0.1324  0.2206    1.0000
3.000000e+07  0.0397  0.0662    0.7086
1.000000e+07  0.0132  0.0221    0.2362
3.000000e+06  0.0040  0.0066    0.0709
1.000000e+06  0.0013  0.0022    0.0236
3.000000e+05  0.0004  0.0007    0.0071
1.000000e+05  0.0001  0.0002    0.0024
```

Each condition reaches baseline at a different budget — richer medium supports more growth but demands more enzyme. Per-condition knee (≥50% baseline): PDA 1.0e+09, half 3.0e+08, onetenth 3.0e+07.

## Budget-limiting enzymes at each condition's knee

Saturated (utilization ≥0.98) — PDA: 0, half: 1, onetenth: 3

Top enzymes by capacity utilization:
```
condition      rxn_id                                                                         name       ec   kcat    flux  utilization
     half rxn09552_c0                  diacylglycerol cholinephosphotransferase, yeast-specific_c0  2.7.8.1  0.065 0.11145       1.0000
 onetenth rxn30674_m0                                                                  Fumarase_m0  4.2.1.2  0.098 0.05922       1.0000
 onetenth PINOS_SC_c0                                Phosphatidylinositol synthase  yeast specific 2.7.8.11  0.065 0.02413       1.0000
 onetenth rxn09524_c0                               NADH dehydrogenase, cytosolic/mitochondrial_m0 3.1.26.5  1.300 0.43939       1.0000
 onetenth rxn03084_c0 5'-Phosphoribosylformylglycinamide:L-glutamine amido-ligase (ADP-forming)_c0  6.3.5.3  0.050 0.00189       0.9592
 onetenth rxn09552_c0                  diacylglycerol cholinephosphotransferase, yeast-specific_c0  2.7.8.1  0.065 0.02386       0.7639
      PDA rxn09552_c0                  diacylglycerol cholinephosphotransferase, yeast-specific_c0  2.7.8.1  0.065 0.33676       0.7550
     half rxn30674_m0                                                                  Fumarase_m0  4.2.1.2  0.098 0.32602       0.7199
 onetenth rxn00533_c0                            Acetyl-CoA:carbon-dioxide ligase (ADP-forming)_c0  6.4.1.2 16.000 0.31472       0.6990
      PDA PINOS_SC_c0                                Phosphatidylinositol synthase  yeast specific 2.7.8.11  0.065 0.34059       0.6886
     half rxn09524_c0                               NADH dehydrogenase, cytosolic/mitochondrial_m0 3.1.26.5  1.300 2.37048       0.5978
      PDA rxn08617_c0              glucose transport via diffusion (extracellular to periplasm)_c0          10.000 5.02656       0.5972
     half PINOS_SC_c0                                Phosphatidylinositol synthase  yeast specific 2.7.8.11  0.065 0.11271       0.5762
      PDA rxn09524_c0                               NADH dehydrogenase, cytosolic/mitochondrial_m0 3.1.26.5  1.300 7.32130       0.5496
 onetenth rxn00216_c0                                                                               2.7.1.1  1.670 0.25089       0.5174
```

## Outputs
- `reaction_enzyme_cost.tsv`, `ec_growth_curve.tsv`, `ec_flux_matrix.tsv`, `ec_saturated_enzymes.tsv`