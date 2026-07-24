#!/opt/env/modelseed/bin/python3
"""Stage 0 — gene → reaction expression aggregation.

Reads the S1 normalized-expression xlsx and the V10 model, aggregates each
reaction's GPR into a single expression score using the standard fungal-GEM
convention:
    score(reaction) = max_over_OR_clauses( min_over_AND_tokens( expr[gene] ) )
i.e. isoenzymes contribute the loudest voice, complex subunits contribute
the quietest.

Emits `outputs/reaction_expression.tsv` (one row per reaction) with:
    rxn_id, name, gpr, n_genes, n_genes_with_expr,
    G1, G2, G3, Mean, agg_mean, agg_G1, agg_G2, agg_G3,
    expression_bin (hi/med/lo/absent), notes.

Bins are set on the Mean log2(TPM+1) distribution of *model genes* (not the
full 13k-gene transcriptome), so bins reflect metabolic-network background:
    hi  ≥ 75th percentile
    med  25th–75th
    lo   > 0 and < 25th
    absent  gene missing from RNA-seq OR aggregate score is nan

Prints a coverage report and fails loudly if the sniff numbers drift.
"""
import json
import re
import sys
from pathlib import Path

import cobra
import numpy as np
import pandas as pd

ROOT = Path('/home/janakae/fungalTemplate/imm904CobraModel/rna_seq_integration')
MODEL_PATH = Path('/home/janakae/fsp237/simulations/gapfill_v1_v2/models/'
                  'fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json')
EXPR_XLSX = ROOT / 'data' / 'S1_normalized_expression.xlsx'
OUT_TSV = ROOT / 'outputs' / 'reaction_expression.tsv'
COVERAGE_JSON = ROOT / 'outputs' / 'coverage_summary.json'

# Expected sniff numbers — verify pipeline reproducibility.
EXPECTED = {
    'n_model_genes': 1274,
    'n_expr_rows': 13047,
    'n_covered': 1174,
    'n_orphan_rxn': 577,
    'n_ch63r_unmigrated': 27,        # gene ids still on the pre-overhaul CH63R_* namespace
    'n_yeast_np_unmigrated': 73,     # legacy yeast (Y*) + RefSeq (NP_*) + SPONT tokens
    'n_model_only_total': 100,       # 27 + 73 (all model genes absent from RNA-seq)
}

TOKEN_RE = re.compile(r'[A-Za-z][A-Za-z0-9_.]*')


def _split_top_level(s, sep):
    parts, depth, start, i = [], 0, 0, 0
    while i < len(s):
        c = s[i]
        if c == '(':
            depth += 1
        elif c == ')':
            depth -= 1
        if depth == 0 and s[i:i + len(sep)] == sep:
            parts.append(s[start:i])
            start = i + len(sep)
            i += len(sep)
            continue
        i += 1
    parts.append(s[start:])
    return parts


def parse_gpr(gpr):
    """Return an OR-of-ANDs list. e.g. '(A and B) or C' → [['A','B'], ['C']]."""
    if not gpr or not gpr.strip():
        return []
    or_clauses = _split_top_level(gpr.strip(), ' or ')
    out = []
    for cl in or_clauses:
        cl = cl.strip()
        if cl.startswith('(') and cl.endswith(')'):
            cl = cl[1:-1].strip()
        toks = [t.strip() for t in re.split(r'\s+and\s+', cl) if t.strip()]
        if toks:
            out.append(toks)
    return out


def aggregate(gpr_str, expr_by_gene):
    """Compute aggregated expression score for a GPR string.

    expr_by_gene: dict gene_id -> numeric score (log2(TPM+1) recommended).

    Returns (score, n_total_genes, n_with_expr). Score is np.nan if the
    reaction has a GPR but no listed gene has an expression measurement.
    Empty GPR returns (np.nan, 0, 0).
    """
    clauses = parse_gpr(gpr_str)
    if not clauses:
        return np.nan, 0, 0
    all_genes = {g for cl in clauses for g in cl}
    n_total = len(all_genes)
    n_with_expr = sum(1 for g in all_genes if g in expr_by_gene)
    if n_with_expr == 0:
        return np.nan, n_total, 0
    or_maxes = []
    for cl in clauses:
        vals = [expr_by_gene[g] for g in cl if g in expr_by_gene]
        if vals:
            or_maxes.append(min(vals))
    if not or_maxes:
        return np.nan, n_total, n_with_expr
    return max(or_maxes), n_total, n_with_expr


def _bin(score, cuts):
    if score is None or (isinstance(score, float) and np.isnan(score)):
        return 'absent'
    if score <= 0:
        return 'absent'
    if score >= cuts['hi']:
        return 'hi'
    if score >= cuts['lo']:
        return 'med'
    return 'lo'


def main():
    model = cobra.io.load_json_model(str(MODEL_PATH))
    model_genes = {g.id for g in model.genes}
    print(f'model : {MODEL_PATH.name}', flush=True)
    print(f'model : {len(model_genes)} genes, {len(model.reactions)} reactions', flush=True)

    xl_log = pd.read_excel(EXPR_XLSX, sheet_name='Normalized log2 TPM')
    xl_tpm = pd.read_excel(EXPR_XLSX, sheet_name='TPM reconstructed')
    xl_log = xl_log.set_index('Gene ID')
    xl_tpm = xl_tpm.set_index('Gene ID')

    expr_rows = set(xl_log.index.astype(str))
    covered = model_genes & expr_rows
    ch63r_unmigrated = {g for g in model_genes if g.startswith('CH63R_')}
    yeast_np_unmigrated = {g for g in (model_genes - expr_rows)
                            if not g.startswith('CH63R_')}
    model_only_total = model_genes - expr_rows

    # Sanity vs sniff
    got = {
        'n_model_genes': len(model_genes),
        'n_expr_rows': len(expr_rows),
        'n_covered': len(covered),
        'n_orphan_rxn': sum(1 for r in model.reactions if not r.gene_reaction_rule.strip()),
        'n_ch63r_unmigrated': len(ch63r_unmigrated),
        'n_yeast_np_unmigrated': len(yeast_np_unmigrated),
        'n_model_only_total': len(model_only_total),
    }
    mismatch = {k: (EXPECTED[k], got[k]) for k in EXPECTED if EXPECTED[k] != got[k]}
    if mismatch:
        print('SNIFF MISMATCH:', mismatch, file=sys.stderr)
        sys.exit(2)
    for k, v in got.items():
        print(f'  {k:22s}: {v}')

    mean_by_gene = xl_log['Mean log2(TPM+1)'].to_dict()
    g1_by_gene = xl_log['G1 log2(TPM+1)'].to_dict()
    g2_by_gene = xl_log['G2 log2(TPM+1)'].to_dict()
    g3_by_gene = xl_log['G3 log2(TPM+1)'].to_dict()
    mean_tpm_by_gene = xl_tpm['Mean TPM'].to_dict()

    # Bin cutoffs on model-gene distribution of Mean log2(TPM+1).
    model_scores = pd.Series({g: mean_by_gene.get(g, np.nan) for g in model_genes}).dropna()
    q25, q75 = float(model_scores.quantile(0.25)), float(model_scores.quantile(0.75))
    cuts = {'hi': q75, 'lo': q25}
    print(f'\nbin cutoffs (Mean log2(TPM+1) over model genes):')
    print(f'  hi  ≥ {q75:.3f}    med [{q25:.3f}, {q75:.3f})    lo (0, {q25:.3f})    absent = 0 / missing')

    rows = []
    for r in model.reactions:
        gpr = r.gene_reaction_rule or ''
        agg_mean, n_tot, n_expr = aggregate(gpr, mean_by_gene)
        agg_g1, _, _ = aggregate(gpr, g1_by_gene)
        agg_g2, _, _ = aggregate(gpr, g2_by_gene)
        agg_g3, _, _ = aggregate(gpr, g3_by_gene)
        agg_tpm, _, _ = aggregate(gpr, mean_tpm_by_gene)
        note = []
        if not gpr.strip():
            note.append('orphan_no_gpr')
        elif n_expr == 0:
            note.append('gpr_present_no_expression')
        if any(t.startswith('CH63R_') for t in TOKEN_RE.findall(gpr)):
            note.append('has_ch63r_unmigrated')
        rows.append({
            'rxn_id': r.id,
            'name': r.name or '',
            'subsystem': r.subsystem or '',
            'gpr': gpr,
            'n_genes': n_tot,
            'n_genes_with_expr': n_expr,
            'agg_mean_log2TPMp1': agg_mean,
            'agg_G1': agg_g1,
            'agg_G2': agg_g2,
            'agg_G3': agg_g3,
            'agg_mean_TPM': agg_tpm,
            'expression_bin': _bin(agg_mean, cuts),
            'notes': ';'.join(note),
        })

    df = pd.DataFrame(rows)
    header = (
        f'# reaction expression scores (Stage 0)\n'
        f'# model = {MODEL_PATH.name}\n'
        f'# expression = {EXPR_XLSX.name} sheet=Normalized log2 TPM\n'
        f'# aggregation = max_over_OR( min_over_AND( gene_expr ) )\n'
        f'# bins on Mean log2(TPM+1) over {len(model_scores)} model genes:'
        f' hi>={q75:.3f}  med=[{q25:.3f},{q75:.3f})  lo=(0,{q25:.3f})  absent=0/missing\n'
    )
    with open(OUT_TSV, 'w') as f:
        f.write(header)
        df.to_csv(f, sep='\t', index=False)
    print(f'\nwrote {OUT_TSV}   ({len(df)} rows)')

    bin_counts = df['expression_bin'].value_counts().to_dict()
    print(f'expression_bin counts: {bin_counts}')

    coverage = {
        'model_genes': len(model_genes),
        'expression_rows': len(expr_rows),
        'model_genes_with_expression': len(covered),
        'ch63r_unmigrated': sorted(ch63r_unmigrated),
        'yeast_np_unmigrated': sorted(yeast_np_unmigrated),
        'orphan_reactions': got['n_orphan_rxn'],
        'reactions_with_gpr': int((df['n_genes'] > 0).sum()),
        'reactions_gpr_no_expression': int(df['notes'].str.contains('gpr_present_no_expression').sum()),
        'bin_cutoffs': cuts,
        'bin_counts': bin_counts,
        'n_model_genes_scored': int(len(model_scores)),
        'percentiles': {
            'p10': float(model_scores.quantile(0.10)),
            'p25': q25,
            'p50': float(model_scores.quantile(0.50)),
            'p75': q75,
            'p90': float(model_scores.quantile(0.90)),
        },
    }
    COVERAGE_JSON.write_text(json.dumps(coverage, indent=2))
    print(f'wrote {COVERAGE_JSON}')


if __name__ == '__main__':
    main()
