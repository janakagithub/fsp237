# Oleate (C18:1, cis-Δ9) β-oxidation in V10

**Model:** `fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json`
**Test condition:** `cpd15269_e0` uptake at 1.0 mmol/gDW/h, aerobic
**Biomass yield:** 0.128 1/h

The V10 model carries TWO routes for oleate β-oxidation in the peroxisome:

1. **A pre-existing LUMPED reaction (rxn09467 / rxn09468)** that converts
   oleoyl-CoA all the way to 9 acetyl-CoA in a single stoichiometry — this
   is what FBA actually uses on the oleate-only test (carries +0.99 flux).
2. **The chain-specific cycle** (C18 → C16 → C14 → … → C4) is also present
   (via V1's gap-fill + the existing C18/C16/C14/C12/C10 cycle), but FBA
   prefers the lumped form because it's the shortest path and both yield
   the same biomass.

Both routes lead to the same net stoichiometry: **1 oleoyl-CoA → 9 acetyl-CoA**
with 8 NADH + 7-8 FADH₂-equivalent (as H₂O₂ via peroxisomal oxidase) +
~7-8 H⁺. The Δ⁹-cis double bond is "consumed" inside the lumped reaction
without needing an explicit Δ³,Δ²-isomerase intermediate.

---

## The path as FBA uses it (lumped, in V10)

```
   ┌────────────────────────── EXTRACELLULAR (e0) ──────────────────────────┐
   │                                                                        │
   │   Oleate                                                               │
   │   cpd15269_e0  ◄── EX_cpd15269_e0  (uptake @ −1.0 mmol/gDW/h)          │
   │                                                                        │
   └────────────────────────────────────────────────────────────────────────┘
                                  │
                                  │  rxn09035_c0
                                  │  Octadecenoate transport
                                  │  (facilitated diffusion, irreversible)
                                  │  octadecenoate_e0 → octadecenoate_c0
                                  │  flux: +1.0
                                  ▼
   ┌────────────────────────────── CYTOSOL (c0) ────────────────────────────┐
   │                                                                        │
   │   Oleate                                                               │
   │   cpd15269_c0                                                          │
   │       │                                                                │
   │       │  rxn08455_c0   fatty-acid-CoA ligase (octadecenoate transport  │
   │       │                via vectoral CoA coupling)                      │
   │       │  ATP + CoA + Oleate ⇌ AMP + PPi + Oleoyl-CoA                   │
   │       │  flux: +1.0  (runs forward = activation)                       │
   │       ▼                                                                │
   │   Oleoyl-CoA                                                           │
   │   cpd15274_c0                                                          │
   │                                                                        │
   └────────────────────────────────────────────────────────────────────────┘
                                  │
                                  │  rxn09844_x0
                                  │  fatty acyl-CoA peroxisomal transport
                                  │  via ABC system (ATP-dependent)
                                  │  H₂O + ATP + Oleoyl-CoA_c0
                                  │  → ADP + Pi + H⁺ + Oleoyl-CoA_x0
                                  │  flux: +0.9945
                                  ▼
   ┌────────────────────────── PEROXISOME (x0) ─────────────────────────────┐
   │                                                                        │
   │   Oleoyl-CoA                                                           │
   │   cpd15274_x0                                                          │
   │       │                                                                │
   │       │  rxn09467_x0   Fatty acid β-oxidation, lumped 8-round          │
   │       │                cycle (oxidase-coupled)                         │
   │       │                                                                │
   │       │      8 H₂O + 8 NAD⁺ + 7 O₂ + 8 CoA + Oleoyl-CoA                │
   │       │      ──────────────────────────────────────────►               │
   │       │      8 NADH + 9 Acetyl-CoA + 7 H₂O₂ + 8 H⁺                     │
   │       │                                                                │
   │       │  (only 7 O₂ instead of 8 because the Δ⁹-cis double bond        │
   │       │   skips one acyl-CoA oxidase step; isomerase auxiliary         │
   │       │   handles Δ³-cis → Δ²-trans conversion that lets the           │
   │       │   shortened chain re-enter the saturated cycle)                │
   │       │                                                                │
   │       │  bounds: (0, ∞) — degradation-direction-only (always was)      │
   │       │  flux: +0.9945                                                 │
   │       │                                                                │
   │       │  Alternative: rxn09468_x0 — same stoichiometry plus 1 NADPH    │
   │       │  consumed by 2,4-dienoyl-CoA reductase auxiliary. FBA picks    │
   │       │  rxn09467 by default (avoids the NADPH cost).                  │
   │       ▼                                                                │
   │   9 Acetyl-CoA_x0                                                      │
   │   (cpd00022_x0)                                                        │
   │       │                                                                │
   │       │  ╭───────────────╮                                             │
   │       │  │ to glyoxylate │   rxn20162_x0  malate synthase              │
   │       │  │   cycle  ──►  │   Acetyl-CoA + Glyoxylate → Malate + CoA    │
   │       │  ╰───────────────╯                                             │
   │       │                                                                │
   └───────┼────────────────────────────────────────────────────────────────┘
           │
           │  (acetyl-CoA also exported to cytosol via existing carnitine
           │   shuttle / direct transport; feeds gluconeogenesis →
           │   biomass precursors)
           ▼
   biomass + H₂O₂ detoxified by peroxisomal catalase
```

### Reactions in this path (in flux order)

| Step | Reaction ID | Compartment | Equation | EC / function | Source |
|---|---|---|---|---|---|
| 0 | `EX_cpd15269_e0` | exchange | Oleate_e0 → ∅ | uptake | existing |
| 1 | `rxn09035_c0` | e0 → c0 | Oleate_e0 → Oleate_c0 | facilitated transport | existing |
| 2 | `rxn08455_c0` | c0 | ATP + CoA + Oleate_c0 ⇌ PPi + AMP + Oleoyl-CoA_c0 | 6.2.1.3 (acyl-CoA ligase) | existing |
| 3 | `rxn09844_x0` | c0 → x0 | H₂O + ATP_x0 + Oleoyl-CoA_c0 → ADP_x0 + Pi_x0 + H⁺_x0 + Oleoyl-CoA_x0 | ABC transporter (Pxa1/Pxa2) | existing |
| 4 | `rxn09467_x0` | x0 | **8 H₂O + 8 NAD + 7 O₂ + 8 CoA + Oleoyl-CoA → 8 NADH + 9 Acetyl-CoA + 7 H₂O₂ + 8 H⁺** | Lumped peroxisomal β-ox of C18:1 (oxidase-coupled, ECI-assisted) | existing |
| 5 | `rxn20162_x0` | x0 | Acetyl-CoA + Glyoxylate → Malate + CoA | 2.3.3.9 (malate synthase) | existing |
| 6 | `rxn00336_c0` | c0 | Isocitrate → Succinate + Glyoxylate | 4.1.3.1 (isocitrate lyase) | existing |

### Alternative oleate-CoA reaction (V10 has both)

| Reaction ID | Equation | Difference |
|---|---|---|
| `rxn09467_x0` | 8 H₂O + 8 NAD + 7 O₂ + 8 CoA + Oleoyl-CoA → 8 NADH + 9 AcCoA + 7 H₂O₂ + 8 H⁺ | "Oxidase route" — uses 2,4-dienoyl path with ECI isomerase only |
| `rxn09468_x0` | 8 H₂O + 8 NAD + **1 NADPH** + 8 O₂ + 8 CoA + Oleoyl-CoA → 8 NADH + 1 NADP + 9 AcCoA + **8** H₂O₂ + 7 H⁺ | "Reductase route" — uses Δ²,⁴-dienoyl-CoA reductase (Sps19) which consumes 1 NADPH; corrects the cis-Δ⁹ via reduction |

In real cells both routes are active and dominance depends on which auxiliary
enzyme (ECI vs DCR) is the rate-limiting step on the particular chain. The
oxidase route (rxn09467) requires only ECI; the reductase route (rxn09468)
covers chains where the Δ³ intermediate becomes Δ²,Δ⁴-dienoyl-CoA after one
more β-ox round and needs the dienoyl-CoA reductase to recover. FBA picks
rxn09467 in V10 because the lumped form skips NADPH — biology is identical
either way (same net 9 acetyl-CoA per oleate).

---

## Chain-specific alternative (also in V10, but not picked by FBA)

The model also has the full chain-by-chain saturated β-oxidation cycle that
could degrade oleoyl-CoA via:

```
Oleoyl-CoA (C18:1) → trans-2,cis-9-Octadecenoyl-CoA → ... 3 rounds of saturated β-ox  →
   cis-3-Dodecenoyl-CoA (C12:1)
        │
        │  Δ³,Δ²-enoyl-CoA isomerase  (ECI)   *** missing as explicit rxn ***
        ▼
   trans-2-Dodecenoyl-CoA (C12:1 → enters saturated cycle)
        │
        │  C12 hydratase + 3-OH-DH + thiolase (existing in V10)
        ▼
   Decanoyl-CoA (C10) → … → Butyryl-CoA (C4) → 2 Acetyl-CoA
```

This chain-specific path would require an explicit Δ³,Δ²-enoyl-CoA
isomerase reaction (cis-3-enoyl → trans-2-enoyl) which is **not currently
in V10**. That's why FBA falls back to the lumped `rxn09467_x0`. If you
ever want to remove the lumped reaction and force chain-specific flux,
you'd need to add `rxn0XXXX_x0` for ECI on each unsaturated intermediate
(C18:1, C16:1, C14:1, C12:1). For now the lumped form is more
parsimonious and produces the same biomass.

---

## Key takeaways

- Oleate is **fully oxidizable aerobically in V10** at 0.128 1/h biomass.
- The actual FBA path is **5 reactions** end-to-end: exchange → transport →
  activation → ABC into peroxisome → lumped β-ox.
- Each oleate yields **9 acetyl-CoA** (vs 8 for palmitate), reflecting the
  one extra carbon pair from C18 vs C16.
- The Δ⁹-cis double bond is handled implicitly inside the lumped reaction;
  the alternative `rxn09468_x0` makes the NADPH cost of the dienoyl-CoA
  reductase explicit but isn't picked by FBA.
- Anaerobically biomass = 0 — β-ox requires O₂ for the acyl-CoA oxidase
  (which makes H₂O₂) and for ETC reoxidation of the 8 NADH produced. No
  anaerobic FA path exists in fungi without acceptor substrates.
