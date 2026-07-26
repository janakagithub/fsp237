#!/opt/env/modelseed/bin/python3
"""Stage 5 — pathway breakdown of the S1 × V10 analyses.

The V10 model only has `reaction.subsystem` set on ~55 gap-filled reactions.
The other 1567 reactions are unassigned, so we impose a heuristic pathway
classification driven by:
  1. reaction.subsystem when set (55 rxns) — highest priority
  2. reaction ID prefix (bio_/EX_/TRP_/tx_ …)
  3. metabolite ID patterns (cpdNNNNN of ATP, AA pool, cofactors, sugars…)
  4. compartment (r0 → ER; x0 → peroxisome)
  5. fallback → "Other / unassigned"

18 pathway buckets — biologically defensible groupings widely used in fungal
GEMs.

Emits four output tables consumed by the site builder:
  outputs/pathway_summary.tsv           one row per pathway
  outputs/pathway_condition_matrix.tsv  pathway × (condition, O2) flux totals
  outputs/reporter_metabolites.tsv      Patil-Nielsen style top metabolites
  outputs/pathway_enrichment.tsv        Fisher's exact for hi-expression enrichment
  outputs/biomass_pathway_corr.tsv      pathway ↔ biomass Pearson r across conditions
  outputs/pathway_assignment.tsv        per-reaction pathway assignment (auditability)
"""
import math
import sys
from collections import defaultdict
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, norm, pearsonr

sys.path.insert(0, '/home/janakae/fsp237/simulations')
from run_simulation_panel import CONDITIONS  # noqa: E402

ROOT = Path('/home/janakae/fungalTemplate/imm904CobraModel/rna_seq_integration')
MODEL_PATH = Path('/home/janakae/fsp237/simulations/gapfill_v1_v2/models/'
                  'fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json')
REXP_TSV = ROOT / 'outputs' / 'reaction_expression.tsv'
STAGE1_DIR = ROOT / 'outputs' / 'stage1_overlay'
STAGE1_SUM = ROOT / 'outputs' / 'stage1_summary.tsv'
OUT_ASSIGN = ROOT / 'outputs' / 'pathway_assignment.tsv'
OUT_SUM = ROOT / 'outputs' / 'pathway_summary.tsv'
OUT_MATRIX = ROOT / 'outputs' / 'pathway_condition_matrix.tsv'
OUT_REPORTER = ROOT / 'outputs' / 'reporter_metabolites.tsv'
OUT_ENRICH = ROOT / 'outputs' / 'pathway_enrichment.tsv'
OUT_BIO_CORR = ROOT / 'outputs' / 'biomass_pathway_corr.tsv'

FLUX_EPS = 1e-6

# --- pathway buckets ------------------------------------------------------
# Order matters — first matching rule wins.

# ModelSEED cpd IDs used by the rules below (stripped of `_<compartment>` suffix)
GLYCOLYSIS_CPDS = {
    'cpd00027',  # glucose
    'cpd00079',  # G6P
    'cpd00089',  # F6P
    'cpd00072',  # F1,6BP  (fructose-1,6-bisphosphate)
    'cpd00102',  # DHAP
    'cpd00169',  # G3P
    'cpd00061',  # 3-PG
    'cpd00169',  # (repeat)
    'cpd00020',  # pyruvate
    'cpd00179',  # maltose (glycogen/starch link)
    'cpd00108',  # galactose
    'cpd00190',  # 1,3-bisphosphoglycerate
    'cpd00069',  # 2-PG (also Tyr) — kept as glycolysis via ID pattern
    'cpd00100',  # glycerol
    'cpd00080',  # PEP
}
PPP_CPDS = {
    'cpd00171',  # 6-P-gluconolactone
    'cpd00072',  # (also FBP)
    'cpd00157',  # 6-phosphogluconate
    'cpd00040',  # ribulose-5-P
    'cpd00108',  # (also galactose)
    'cpd00284',  # xylulose-5-P
    'cpd00280',  # (also galU) — do NOT include here; pattern-guarded below
    'cpd00214',  # (palmitate — exclude)
    'cpd00212',  # ribose-5-P
    'cpd00227',  # sedoheptulose-7-P
    'cpd00236',  # erythrose-4-P
}
TCA_CPDS = {
    'cpd00137',  # citrate
    'cpd00417',  # cis-aconitate
    'cpd00157',  # (also 6PG — pattern-guarded)
    'cpd00417',  # (repeat)
    'cpd00024',  # 2-oxoglutarate (α-KG)
    'cpd00036',  # succinate
    'cpd00106',  # fumarate
    'cpd00130',  # malate
    'cpd00032',  # oxaloacetate
    'cpd00081',  # succinyl-CoA
}
AA_20 = {
    'cpd00023', 'cpd00033', 'cpd00035', 'cpd00039', 'cpd00041', 'cpd00051',
    'cpd00053', 'cpd00054', 'cpd00060', 'cpd00065', 'cpd00066', 'cpd00069',
    'cpd00084', 'cpd00107', 'cpd00119', 'cpd00129', 'cpd00132', 'cpd00156',
    'cpd00161', 'cpd00322',
}
NUC_TRIPHOS = {'cpd00002', 'cpd00038', 'cpd00052', 'cpd00062'}  # ATP, GTP, CTP, UTP
NUC_DIPHOS  = {'cpd00008', 'cpd00031', 'cpd00046', 'cpd00093'}  # ADP, GDP, CDP, UDP
COFACTORS   = {
    'cpd00003', 'cpd00004',   # NAD+, NADH
    'cpd00005', 'cpd00006',   # NADPH, NADP+
    'cpd00010',               # CoA
    'cpd00015',               # FAD
    'cpd00982',               # FADH2
    'cpd00025',               # H2O2
    'cpd00013',               # NH3
    'cpd00007',               # O2
}
CELL_WALL = {
    'cpd11683',  # chitin
    'cpd11791',  # β-1,3-glucan
    'cpd12148',  # α-1,3-glucan
    'cpd11685',  # mannan
    'cpd00135',  # UDP-GlcNAc
    'cpd00232',  # UDP-glucose
}
STORAGE = {
    'cpd00155',  # glycogen
    'cpd00794',  # trehalose
}
PECTIN = {
    'cpd00280',  # D-galacturonate
    'cpd11454',  # L-galactonate
    'cpd11458',  # (2-dehydro-3-deoxy-L-galactonate)
}
PENTOSE_SUGAR = {
    'cpd00154',  # xylose
    'cpd00224',  # L-arabinose
    'cpd00379',  # L-arabitol
    'cpd00340',  # L-xylulose
    'cpd00306',  # xylitol
}
FA_LONG = {
    'cpd00214',  # palmitate
    'cpd15269',  # oleate
    'cpd15240',  # hexacosanoate
}
FA_LONG_TRIGGER_SUFFIX = 'CoA'  # heuristic: reaction name mentions CoA + long chain

# ---------- classification --------------------------------------------------

def _strip_comp(mid: str) -> str:
    return mid.rsplit('_', 1)[0] if '_' in mid else mid


def classify_pathway(rxn) -> str:
    """Return the pathway bucket a reaction belongs to.
    Rules cascade top-down — first hit wins."""

    rid = rxn.id
    subs = (rxn.subsystem or '').strip().lower()
    met_ids = [_strip_comp(m.id) for m in rxn.metabolites]
    met_set = set(met_ids)
    compartments = {m.id.rsplit('_', 1)[-1] for m in rxn.metabolites if '_' in m.id}
    name_lower = (rxn.name or '').lower()

    # 0. explicit
    if subs:
        if 'beta-oxidation' in subs or subs.startswith('beta-ox'):
            return 'β-oxidation (peroxisomal)'
        if 'peroxisomal cofactor' in subs:
            return 'β-oxidation (peroxisomal)'
        if 'vlcfa' in subs:
            return 'β-oxidation (peroxisomal)'
        if 'ashwell' in subs or 'galu' in subs:
            return 'Pectin / D-galacturonate (Ashwell)'
        if 'penttilae' in subs or 'arabinose' in subs:
            return 'Pentose-sugar catabolism'
        if 'chitin' in subs or 'alpha-1,3-glucan' in subs or 'cell wall' in subs:
            return 'Cell-wall polysaccharide'
        if 'mannitol' in subs:
            return 'Compatible solute / mannitol'
        if 'melanin' in subs:
            return 'Melanin / DHN'
        # Fall through if subs was set but generic

    # 1. biomass & exchange
    if rid.startswith('bio') or 'biomass' in name_lower:
        return 'Biomass'
    if rid.startswith('EX_'):
        return 'Exchange'
    if rid.startswith('DM_') or rid.startswith('SK_'):
        return 'Sink / demand'

    # 2. Transport
    if rid.startswith('TRP_') or rid.startswith('tx_') or 'transporter' in name_lower or 'transport' in name_lower:
        return 'Transport (inter-compartmental)'

    # 3. Pentose sugar catabolism (before generic glycolysis to catch xylose/arabinose)
    if met_set & PENTOSE_SUGAR:
        return 'Pentose-sugar catabolism'

    # 4. Pectin / Ashwell
    if met_set & PECTIN:
        return 'Pectin / D-galacturonate (Ashwell)'

    # 5. Storage carbohydrates
    if met_set & STORAGE:
        return 'Storage (glycogen / trehalose)'

    # 6. Cell wall
    if met_set & CELL_WALL:
        return 'Cell-wall polysaccharide'

    # 7. β-oxidation — peroxisomal + fatty-acyl-CoA or long FA name
    if 'x0' in compartments and (any(k in name_lower for k in ('acyl-coa', 'ketoacyl', 'enoyl', 'hydroxyacyl'))
                                 or met_set & FA_LONG):
        return 'β-oxidation (peroxisomal)'

    # 8. Fatty-acid biosynthesis (cytosolic malonyl/acyl-ACP)
    if 'malonyl' in name_lower or 'acp' in name_lower or 'acyl carrier' in name_lower:
        return 'Fatty-acid biosynthesis'
    if met_set & FA_LONG and 'x0' not in compartments:
        return 'Fatty-acid biosynthesis'

    # 9. TCA
    if len(met_set & TCA_CPDS) >= 2 or ('m0' in compartments and met_set & TCA_CPDS):
        # need at least one specifically TCA-diagnostic met
        if met_set & {'cpd00137', 'cpd00417', 'cpd00024', 'cpd00036',
                       'cpd00106', 'cpd00130', 'cpd00032', 'cpd00081'}:
            return 'TCA cycle'

    # 10. ETC / OxPhos — reactions moving H+ between compartments + cofactor
    if ('cpd00067' in met_set) and (len({'m0', 'c0'} & compartments) == 2):
        return 'Oxidative phosphorylation / ETC'
    if ('cpd00007' in met_set and 'cpd00067' in met_set) and 'm0' in compartments:
        return 'Oxidative phosphorylation / ETC'

    # 11. Amino acid metabolism
    if met_set & AA_20:
        return 'Amino-acid metabolism'

    # 12. Nucleotide metabolism (must have NTP/NDP as substrate AND a
    # nucleobase-like companion — otherwise ATP hydrolysis dominates)
    if (met_set & NUC_TRIPHOS) and any(cpd in met_set for cpd in (
        'cpd00046', 'cpd00126', 'cpd00062', 'cpd00091',  # UMP/GMP/UTP/dGMP
        'cpd00206', 'cpd00294', 'cpd00296', 'cpd00298',
    )):
        return 'Nucleotide metabolism'

    # 13. PPP — key diagnostic mets
    if met_set & {'cpd00171', 'cpd00157', 'cpd00040', 'cpd00284',
                   'cpd00212', 'cpd00227', 'cpd00236'}:
        return 'Pentose phosphate pathway'

    # 14. Glycolysis / gluconeogenesis
    if met_set & GLYCOLYSIS_CPDS:
        return 'Glycolysis / gluconeogenesis'

    # 15. Cofactor / vitamin biosynthesis
    if any(k in name_lower for k in ('nad', 'coenzyme a', 'fmn', 'thiamin', 'folate',
                                       'biotin', 'riboflavin', 'pantothenate', 'pyridoxal')):
        return 'Cofactor / vitamin biosynthesis'

    # 16. Sulfur / nitrogen
    if any(cpd in met_set for cpd in ('cpd00048', 'cpd00013')) and (
       'sulf' in name_lower or 'ammon' in name_lower or 'nitrat' in name_lower or 'urea' in name_lower
    ):
        return 'Sulfur / nitrogen assimilation'

    # 17. Lipid / membrane (glycerophospholipid, sterol)
    if any(k in name_lower for k in ('phosphatidyl', 'phosphocholine', 'ergosterol', 'sphingo',
                                       'cardiolipin', 'lanosterol', 'zymosterol')):
        return 'Lipid / membrane'

    # 18. GAM stoich (ATP hydrolysis w/ H2O)
    if met_set == {'cpd00001', 'cpd00002', 'cpd00008', 'cpd00009', 'cpd00067'}:
        return 'GAM / maintenance'

    return 'Other / unassigned'


# ---------- reporter metabolites (Patil-Nielsen) --------------------------

def reporter_metabolites(model, rxn_score):
    """For each metabolite, aggregate scores of its neighbouring reactions:
        z_met = mean(z_r for r in reactions_of_met) * sqrt(k_met)
    where z_r = (rxn_score_r - μ) / σ over all scored reactions.

    Returns a DataFrame sorted by |z_met| desc.
    """
    scores = np.array([s for s in rxn_score.values()
                        if s is not None and not math.isnan(s)])
    if len(scores) < 20:
        return pd.DataFrame()
    mu = float(scores.mean())
    sigma = float(scores.std(ddof=1) or 1.0)

    rows = []
    for met in model.metabolites:
        vals = []
        for r in met.reactions:
            s = rxn_score.get(r.id)
            if s is None or (isinstance(s, float) and math.isnan(s)):
                continue
            vals.append((s - mu) / sigma)
        if len(vals) < 3:
            continue
        z_met = float(np.mean(vals) * math.sqrt(len(vals)))
        rows.append({
            'metabolite_id': met.id,
            'name': met.name or '',
            'compartment': met.id.rsplit('_', 1)[-1] if '_' in met.id else '',
            'formula': met.formula or '',
            'k_reactions': len(vals),
            'z_reporter': round(z_met, 3),
            'p_one_sided': round(1 - float(norm.cdf(z_met)), 4) if z_met > 0 else
                            round(float(norm.cdf(z_met)), 4),
        })
    return pd.DataFrame(rows).sort_values('z_reporter', key=lambda c: c.abs(),
                                             ascending=False)


# ---------- main ----------------------------------------------------------

def main():
    print(f'model : {MODEL_PATH.name}')
    model = cobra.io.load_json_model(str(MODEL_PATH))
    rexp = pd.read_csv(REXP_TSV, sep='\t', comment='#').set_index('rxn_id')

    # 1. Assign pathway to every reaction
    assignment = []
    for r in model.reactions:
        pw = classify_pathway(r)
        assignment.append({
            'rxn_id': r.id, 'name': r.name or '', 'subsystem': r.subsystem or '',
            'pathway': pw,
            'compartment_bucket': ','.join(sorted({m.id.rsplit('_',1)[-1] for m in r.metabolites if '_' in m.id})),
        })
    assign_df = pd.DataFrame(assignment)
    assign_df.to_csv(OUT_ASSIGN, sep='\t', index=False)
    print(f'wrote {OUT_ASSIGN.name}')

    n_per_pathway = assign_df.groupby('pathway').size().sort_values(ascending=False)
    print('\npathway assignment counts:')
    for p, n in n_per_pathway.items():
        print(f'  {n:4d}  {p}')

    # 2. Per-pathway expression stats
    merged = assign_df.set_index('rxn_id').join(
        rexp[['agg_mean_log2TPMp1', 'expression_bin', 'n_genes', 'n_genes_with_expr']]
    )
    merged['agg_mean_log2TPMp1'] = pd.to_numeric(merged['agg_mean_log2TPMp1'], errors='coerce')

    pathway_stats = []
    for pw, g in merged.groupby('pathway'):
        scored = g['agg_mean_log2TPMp1'].dropna()
        bin_counts = g['expression_bin'].value_counts().to_dict()
        pathway_stats.append({
            'pathway': pw,
            'n_reactions': int(len(g)),
            'n_with_gpr': int((g['n_genes'] > 0).sum()),
            'n_scored': int(len(scored)),
            'mean_expression': round(float(scored.mean()), 3) if len(scored) else None,
            'median_expression': round(float(scored.median()), 3) if len(scored) else None,
            'hi_n': int(bin_counts.get('hi', 0)),
            'med_n': int(bin_counts.get('med', 0)),
            'lo_n': int(bin_counts.get('lo', 0)),
            'absent_n': int(bin_counts.get('absent', 0)),
            'hi_frac': round(bin_counts.get('hi', 0) / max(1, len(g)), 3),
        })
    stats_df = pd.DataFrame(pathway_stats).sort_values('mean_expression', ascending=False)

    # 3. Fisher enrichment for hi-expression per pathway (vs rest of model)
    total = len(merged)
    total_hi = int((merged['expression_bin'] == 'hi').sum())
    enrich_rows = []
    for pw, g in merged.groupby('pathway'):
        pw_n = int(len(g))
        pw_hi = int((g['expression_bin'] == 'hi').sum())
        other_n = total - pw_n
        other_hi = total_hi - pw_hi
        table = [[pw_hi, pw_n - pw_hi], [other_hi, other_n - other_hi]]
        try:
            odds_r, p = fisher_exact(table, alternative='greater')
        except Exception:
            odds_r, p = float('nan'), 1.0
        enrich_rows.append({
            'pathway': pw, 'n_reactions': pw_n, 'n_hi_in_pathway': pw_hi,
            'expected_hi': round(pw_n * total_hi / max(1, total), 2),
            'odds_ratio': round(float(odds_r) if odds_r == odds_r else 0, 3),
            'p_fisher_greater': round(float(p), 5),
        })
    enrich_df = pd.DataFrame(enrich_rows).sort_values('p_fisher_greater')
    # BH correction
    m_tests = len(enrich_df)
    enrich_df['bh_q'] = 1.0
    if m_tests > 0:
        ps = enrich_df['p_fisher_greater'].values
        ranks = np.arange(1, m_tests + 1)
        bh = ps * m_tests / ranks
        bh = np.minimum.accumulate(bh[::-1])[::-1]
        enrich_df['bh_q'] = np.round(np.clip(bh, 0, 1), 4)
    # Merge enrichment into stats output
    stats_df = stats_df.merge(enrich_df[['pathway', 'odds_ratio', 'p_fisher_greater', 'bh_q']],
                                on='pathway')
    stats_df.to_csv(OUT_SUM, sep='\t', index=False)
    enrich_df.to_csv(OUT_ENRICH, sep='\t', index=False)
    print(f'wrote {OUT_SUM.name}')
    print(f'wrote {OUT_ENRICH.name}')
    print(f'\ntop-8 enriched pathways for hi-expression:')
    for _, r in enrich_df.head(8).iterrows():
        print(f'  q={r["bh_q"]:.4f}  OR={r["odds_ratio"]:>5.2f}  '
              f'{int(r["n_hi_in_pathway"])}/{int(r["n_reactions"])}  {r["pathway"]}')

    # 4. Pathway × condition activity matrix (from Stage 1 pFBA outputs)
    cond_rows = []
    pathway_flux = defaultdict(dict)
    pathway_active = defaultdict(dict)
    for cond_id, label, stage, cs, notes in CONDITIONS:
        for aerobic in (True, False):
            o2 = 'aerobic' if aerobic else 'anaerobic'
            path = STAGE1_DIR / f'{cond_id}_{o2}.tsv'
            if not path.exists():
                continue
            df = pd.read_csv(path, sep='\t')
            df = df.set_index('rxn_id').join(assign_df.set_index('rxn_id')[['pathway']])
            for pw, g in df.groupby('pathway'):
                flux_sum = float(g['flux_pFBA'].abs().sum())
                n_active = int((g['flux_pFBA'].abs() > FLUX_EPS).sum())
                pathway_flux[pw][f'{cond_id}_{o2[:3]}'] = round(flux_sum, 3)
                pathway_active[pw][f'{cond_id}_{o2[:3]}'] = n_active
                cond_rows.append({
                    'pathway': pw, 'condition_id': cond_id, 'stage': stage,
                    'O2': o2, 'sum_abs_flux': round(flux_sum, 3),
                    'n_active': n_active, 'n_total': int(len(g)),
                    'active_rate': round(n_active / max(1, len(g)), 3),
                })
    cond_df = pd.DataFrame(cond_rows)
    cond_df.to_csv(OUT_MATRIX, sep='\t', index=False)
    print(f'wrote {OUT_MATRIX.name}')

    # 5. Biomass ↔ pathway flux correlation across the 18 aerobic conditions
    st1_summary = pd.read_csv(STAGE1_SUM, sep='\t')
    aer = st1_summary[st1_summary['O2'] == 'aerobic'][['condition_id', 'biomass']]
    corr_rows = []
    for pw, g in cond_df[cond_df['O2'] == 'aerobic'].groupby('pathway'):
        joined = g.merge(aer, on='condition_id')
        if joined['sum_abs_flux'].std(ddof=1) == 0 or joined['biomass'].std(ddof=1) == 0:
            continue
        r, p = pearsonr(joined['sum_abs_flux'], joined['biomass'])
        corr_rows.append({
            'pathway': pw, 'n_conditions': len(joined),
            'pearson_r_flux_vs_biomass': round(float(r), 3),
            'p_value': round(float(p), 4),
        })
    corr_df = pd.DataFrame(corr_rows).sort_values(
        'pearson_r_flux_vs_biomass', ascending=False
    )
    corr_df.to_csv(OUT_BIO_CORR, sep='\t', index=False)
    print(f'wrote {OUT_BIO_CORR.name}')
    print(f'\ntop-5 pathways positively correlating with biomass (aerobic):')
    for _, r in corr_df.head(5).iterrows():
        print(f'  r={r["pearson_r_flux_vs_biomass"]:+.3f}  '
              f'p={r["p_value"]:.4f}  {r["pathway"]}')

    # 6. Reporter metabolites (Patil-Nielsen) — S1 expression
    rxn_score = rexp['agg_mean_log2TPMp1'].to_dict()
    rep = reporter_metabolites(model, rxn_score)
    rep.to_csv(OUT_REPORTER, sep='\t', index=False)
    print(f'wrote {OUT_REPORTER.name}   ({len(rep)} metabolites scored)')
    print(f'\ntop-8 reporter metabolites (|z| desc):')
    for _, r in rep.head(8).iterrows():
        print(f'  z={r["z_reporter"]:+.2f}  k={int(r["k_reactions"])}  '
              f'{r["metabolite_id"]:20s}  {r["name"][:60]}')


if __name__ == '__main__':
    main()
