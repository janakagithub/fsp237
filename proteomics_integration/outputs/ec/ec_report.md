# Phase 2 — proteome-allocated enzyme-capacity model

Model fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json. Catalyzed reactions costed: **1045**. kcat = EC-class prior (Bar-Even 2011; median 10.0 s⁻¹, refined by EC first digit) — order-of-magnitude, swappable for BRENDA/DLKcat. MW from FSP237 proteome sequences (median 44.1 kDa).

## Enzyme-budget growth curve (frac of own-medium baseline)

```
condition        PDA    half  onetenth
budget                                
1.000000e+10  1.0000  1.0000    1.0000
3.000000e+09  1.0000  1.0000    1.0000
1.000000e+09  1.0000  1.0000    1.0000
3.000000e+08  0.5240  1.0000    1.0000
1.000000e+08  0.2014  0.6752    1.0000
3.000000e+07  0.0832  0.2312    1.0000
1.000000e+07  0.0484  0.0995    0.5913
3.000000e+06  0.0167  0.0416    0.1774
1.000000e+06  0.0071  0.0154    0.0591
3.000000e+05  0.0036  0.0062    0.0177
1.000000e+05  0.0012  0.0027    0.0059
```

Each condition reaches baseline at a different budget — richer medium supports more growth but demands more enzyme. Per-condition knee (≥50% baseline): PDA 3.0e+08, half 1.0e+08, onetenth 1.0e+07.

## Budget-limiting enzymes at each condition's knee

Saturated (utilization ≥0.98) — PDA: 1, half: 1, onetenth: 1

Top enzymes by capacity utilization:
```
condition      rxn_id                                                            name       ec  kcat     flux  utilization
      PDA rxn08617_c0 glucose transport via diffusion (extracellular to periplasm)_c0          10.00  2.52498       1.0000
     half rxn08617_c0 glucose transport via diffusion (extracellular to periplasm)_c0          10.00  1.65972       1.0000
 onetenth rxn00533_c0               Acetyl-CoA:carbon-dioxide ligase (ADP-forming)_c0 2.1.3.15 28.00  0.26265       1.0000
 onetenth rxn08617_c0 glucose transport via diffusion (extracellular to periplasm)_c0          10.00  0.28884       0.8198
 onetenth rxn09780_c0                           ADP/ATP transporter, mitochondrial_c0          10.00  1.38968       0.5390
     half rxn09780_c0                           ADP/ATP transporter, mitochondrial_c0          10.00  7.94345       0.4457
 onetenth rxn01106_c0                      2-Phospho-D-glycerate 2,3-phosphomutase_c0 5.4.2.12 12.00 -0.38230       0.4150
      PDA rxn09780_c0                           ADP/ATP transporter, mitochondrial_c0          10.00 12.30208       0.2955
 onetenth rxn00168_c0                                                                  1.8.3.2 23.57  0.29130       0.2825
      PDA rxn00533_c0               Acetyl-CoA:carbon-dioxide ligase (ADP-forming)_c0 2.1.3.15 28.00  2.32730       0.1464
 onetenth rxn00175_c0                             Acetate:CoA ligase (AMP-forming)_c0 6.2.1.17 15.00  0.29142       0.1444
     half rxn00533_c0               Acetyl-CoA:carbon-dioxide ligase (ADP-forming)_c0 2.1.3.15 28.00  1.49939       0.1396
 onetenth rxn09524_c0                  NADH dehydrogenase, cytosolic/mitochondrial_m0  1.6.5.9 27.90  0.42004       0.1336
 onetenth rxn00216_c0                                                                  2.7.1.1 28.00  0.28884       0.1066
     half rxn01106_c0                      2-Phospho-D-glycerate 2,3-phosphomutase_c0 5.4.2.12 12.00 -2.20381       0.0999
```

## Outputs
- `reaction_enzyme_cost.tsv`, `ec_growth_curve.tsv`, `ec_flux_matrix.tsv`, `ec_saturated_enzymes.tsv`