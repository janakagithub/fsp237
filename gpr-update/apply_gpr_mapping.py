#!/usr/bin/env python3
"""Update FSP237 model GPRs: substitute CH63R_* tokens with mapped gene_*
identifiers from BLAST results.

Inputs (in /gpr-update/):
  - ch63r_to_fsp237_mapping.tsv   (output of BLAST analysis)
  - ../fsp237_biomass_extension/fsp237_atp_safe_gsm_extended.json  (input model)

Outputs (in /gpr-update/):
  - fsp237_atp_safe_gsm_gpr_updated.json   (updated model)
  - gpr_change_log.tsv                      (per-reaction before/after)
  - flagged_reactions.tsv                   (reactions needing manual review)

Rules:
  - Only touch GPRs that contain at least one CH63R_* token
  - Existing gene_* tokens are preserved
  - Each CH63R_* is replaced by its mapped gene_* (if mapping passes thresholds)
  - Unmapped CH63R_* is REMOVED and the reaction is flagged for review
  - Duplicate gene_* in an "or" clause are deduped (no "gene_X or gene_X")
  - Preserves "or" / "and" logical structure
"""
import json
import os
import re
import sys

import cobra
import pandas as pd

GPR_DIR = '/home/janakae/fungalTemplate/imm904CobraModel/gpr-update'
INPUT_MODEL  = '/home/janakae/fungalTemplate/imm904CobraModel/fsp237_biomass_extension/fsp237_atp_safe_gsm_extended.json'
MAPPING_TSV  = f'{GPR_DIR}/ch63r_to_fsp237_mapping.tsv'
OUTPUT_MODEL = f'{GPR_DIR}/fsp237_atp_safe_gsm_gpr_updated.json'
CHANGELOG    = f'{GPR_DIR}/gpr_change_log.tsv'
FLAGGED      = f'{GPR_DIR}/flagged_reactions.tsv'

# Token regex matches CH63R_NNNN, gene_NNNN (or any other word-like identifier)
TOKEN_RE = re.compile(r'[A-Za-z][A-Za-z0-9_.]*')


def load_mapping():
    df = pd.read_csv(MAPPING_TSV, sep='\t')
    # CH63R -> gene_* (only where we have a mapping)
    m = {row['ch63r_id']: row['mapped_gene']
         for _, row in df.iterrows()
         if isinstance(row['mapped_gene'], str) and row['mapped_gene'].startswith('gene_')}
    weak = {row['ch63r_id'] for _, row in df.iterrows()
            if isinstance(row.get('flag', ''), str) and 'WEAK' in row['flag']}
    unmapped = {row['ch63r_id'] for _, row in df.iterrows()
                if not (isinstance(row['mapped_gene'], str) and row['mapped_gene'].startswith('gene_'))}
    return m, weak, unmapped


def parse_gpr(gpr):
    """Return list of and-clauses (each is a list of OR'd tokens).
    A GPR like '(A and B) or C or (D and E)' becomes [['A','B'], ['C'], ['D','E']].
    A GPR like 'A or B or C' becomes [['A'], ['B'], ['C']].
    A GPR like 'A and B' becomes [['A','B']].
    """
    if not gpr or not gpr.strip():
        return []
    # Normalize whitespace and parentheses spacing
    s = gpr.strip()
    # Split top-level 'or' (respecting parentheses)
    or_clauses = _split_top_level(s, ' or ')
    out = []
    for cl in or_clauses:
        cl = cl.strip()
        if cl.startswith('(') and cl.endswith(')'):
            cl = cl[1:-1].strip()
        and_tokens = [t.strip() for t in re.split(r'\s+and\s+', cl)]
        # filter tokens
        and_tokens = [t for t in and_tokens if t]
        if and_tokens: out.append(and_tokens)
    return out


def _split_top_level(s, sep):
    """Split s on top-level occurrences of sep (not inside parentheses)."""
    parts = []
    depth = 0
    start = 0
    i = 0
    while i < len(s):
        c = s[i]
        if c == '(': depth += 1
        elif c == ')': depth -= 1
        if depth == 0 and s[i:i+len(sep)] == sep:
            parts.append(s[start:i])
            start = i + len(sep)
            i += len(sep)
            continue
        i += 1
    parts.append(s[start:])
    return parts


def build_gpr(clauses):
    """Inverse of parse_gpr."""
    if not clauses:
        return ''
    rendered = []
    for cl in clauses:
        # Dedupe within an 'and' clause
        cl = list(dict.fromkeys(cl))
        if len(cl) == 1:
            rendered.append(cl[0])
        else:
            rendered.append('(' + ' and '.join(cl) + ')')
    # Dedupe across 'or' clauses (each unique whole clause-string)
    rendered = list(dict.fromkeys(rendered))
    if len(rendered) == 1:
        # If the only clause is a parenthesized 'and', strip parens
        only = rendered[0]
        if only.startswith('(') and only.endswith(')'):
            return only[1:-1]
        return only
    return ' or '.join(rendered)


def remap_gpr(gpr, mapping):
    """Apply CH63R->gene_* mapping to a GPR.
    Returns (new_gpr, change_info dict).
    change_info: removed_ch63r=[], replaced=[(ch63r, gene), ...], unmapped_present=[]
    """
    clauses = parse_gpr(gpr)
    replaced, removed = [], []
    new_clauses = []
    for cl in clauses:
        new_tokens = []
        for tok in cl:
            if tok.startswith('CH63R_'):
                if tok in mapping:
                    new_tokens.append(mapping[tok])
                    replaced.append((tok, mapping[tok]))
                else:
                    removed.append(tok)  # unmapped: drop
            else:
                new_tokens.append(tok)
        if new_tokens:
            new_clauses.append(new_tokens)
    new_gpr = build_gpr(new_clauses)
    return new_gpr, {'replaced': replaced, 'removed': removed}


def main():
    print(f'loading: {INPUT_MODEL}')
    m = cobra.io.load_json_model(INPUT_MODEL)
    print(f'  reactions: {len(m.reactions)}, genes: {len(m.genes)}')

    mapping, weak, unmapped = load_mapping()
    print(f'\nmapping table: {len(mapping)} CH63R->gene_* mappings')
    print(f'  weak (pident<50): {sorted(weak)}')
    print(f'  unmapped CH63R  : {sorted(unmapped)}')

    print(f'\n=== UPDATING GPRs ===')
    changelog = []
    flagged = []
    for r in m.reactions:
        gpr = r.gene_reaction_rule or ''
        if 'CH63R_' not in gpr:
            continue
        new_gpr, info = remap_gpr(gpr, mapping)
        old_gpr = gpr
        # Apply
        r.gene_reaction_rule = new_gpr
        changelog.append({
            'rxn_id': r.id, 'name': r.name,
            'old_gpr': old_gpr, 'new_gpr': new_gpr,
            'n_replaced': len(info['replaced']),
            'n_removed_unmapped': len(info['removed']),
            'replaced': ';'.join(f'{a}->{b}' for a,b in info['replaced']),
            'removed': ';'.join(info['removed']),
        })
        flag_msg = []
        if info['removed']:
            flag_msg.append(f'unmapped CH63R removed: {info["removed"]}')
        weak_in = [t for t,_ in info['replaced'] if t in weak]
        if weak_in:
            flag_msg.append(f'weak BLAST matches used: {weak_in}')
        if not new_gpr:
            flag_msg.append('GPR became EMPTY (no genes survive)')
        if flag_msg:
            flagged.append({
                'rxn_id': r.id, 'name': r.name,
                'old_gpr': old_gpr, 'new_gpr': new_gpr,
                'issue': ' | '.join(flag_msg),
            })
        print(f'  {r.id}: {len(info["replaced"])} replaced, {len(info["removed"])} removed (unmapped)')

    # Save model
    cobra.io.save_json_model(m, OUTPUT_MODEL)
    print(f'\nwrote updated model: {OUTPUT_MODEL}')

    # Save change log
    pd.DataFrame(changelog).to_csv(CHANGELOG, sep='\t', index=False)
    print(f'wrote {len(changelog)} rxn changes to: {CHANGELOG}')

    pd.DataFrame(flagged).to_csv(FLAGGED, sep='\t', index=False)
    print(f'wrote {len(flagged)} flagged reactions to: {FLAGGED}')

    # Summary
    print('\n=== SUMMARY ===')
    print(f'  reactions updated     : {len(changelog)}')
    print(f'  reactions flagged     : {len(flagged)}')
    print(f'  CH63R replacements    : {sum(c["n_replaced"] for c in changelog)}')
    print(f'  CH63R unmapped (removed): {sum(c["n_removed_unmapped"] for c in changelog)}')

    # Validate: no CH63R left in any GPR
    leftover = [r.id for r in m.reactions if 'CH63R_' in (r.gene_reaction_rule or '')]
    print(f'  reactions still containing CH63R: {len(leftover)}')
    if leftover:
        for rid in leftover: print(f'    {rid}')

    # Re-run a smoke FBA to ensure model still works
    print('\n=== model smoke-test ===')
    sol = m.optimize()
    print(f'  optimize() returned: {sol.objective_value}')

if __name__ == '__main__':
    main()
