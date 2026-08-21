# PDB Media — Concentration→Uptake Scaling Audit

## Anchor and formula
- **Glucose (dextrose)** is the anchor carbon source. Baseline max glucose uptake = **5.0 mmol gDW⁻¹ h⁻¹**, at PDB glucose concentration **111.92 mM**.
- Base scaling for every nutrient *i*: `uptake_i = glucose_flux × (conc_i / conc_glucose)` = `5.0 × (conc_i / 111.92)`.

## Secondary-carbon cap (κ)
Amino acids and the secondary sugars (fructose, sucrose) contain carbon and, left uncapped,
would supply ~10% of organic carbon. To keep glucose **overwhelmingly dominant** while
**preserving the measured relative abundance** among the secondary sources, a single cap
factor κ is applied to all secondary *organic carbon* sources:

- Uncapped aggregate secondary organic carbon = **3.365 mmol C** (vs glucose 30.0 mmol C).
- κ chosen so aggregate secondary organic C ≤ **5% of glucose carbon** (0.05 × 30 = 1.5 mmol C):
  **κ = min(1, 1.5 / 3.365) = 0.445797**.
- Applied as `uptake_i = 5.0 × (conc_i/111.92) × κ` for fructose, sucrose, and all 20 amino acids.
- Because κ is a single proportional factor, the *relative* abundances measured in PDB are
  preserved exactly; only the absolute secondary-carbon budget is capped.

## Vitamins (trace)
B1, B2, B3, B5, B6 (pyridoxal), B9 are concentration-scaled **without** κ (their carbon is
negligible). Resulting bounds range ~3×10⁻⁶ to 7×10⁻⁴ mmol gDW⁻¹ h⁻¹ — genuinely trace, and
per-compound (not a blanket flat value), driven by each vitamin's own PDB concentration.

## Minerals
- **Phosphate (cpd00009), K⁺ (cpd00205), Na⁺ (cpd00971), Fe²⁺ (cpd10515)** map to the model's
  standard inorganic minimal-media backbone, which every one of the 18 existing panel conditions
  opens at −1000. They are **not** concentration-capped: they are non-carbon, biomass-essential,
  and capping them would confound the carbon-driven dilution response (growth would become
  mineral-limited rather than carbon-limited). Their PDB presence is recorded for provenance in
  the mapping table; the appropriate **phosphate species** is used for phosphorus, not elemental P.
- **Ca, Mg, Zn, Cu, Mn, Vitamin C**: no exchange / no consumption pathway in the model — **omitted
  entirely** (per instruction: compounds the model cannot consume are left out rather than
  force-mapped or gap-filled).

## Dilution series
The three media are the baseline bounds scaled by **1.0 / 0.5 / 0.1** for all PDB-supplied
organics and vitamins (glucose 5.0 / 2.5 / 0.5). The inorganic backbone is held constant at
−1000 across all three (it represents buffered salt capacity, and is constant in every panel
condition).

## Validation (pFBA, CPLEX, objective = bio_gsm)
| Medium | biomass | glucose uptake | glucose C fraction |
|---|---|---|---|
| PDB-baseline | 0.19147 | 5.000 | 95.23% |
| PDB-half | 0.09574 | 2.500 | 95.23% |
| PDB-onetenth | 0.01915 | 0.500 | 95.23% |

Dilution response strictly monotonic (baseline > half > onetenth); glucose is the dominant
carbon source in all three; **no biomass coefficients, ATP stoichiometry, or central-carbon /
ETC reactions were modified**.
