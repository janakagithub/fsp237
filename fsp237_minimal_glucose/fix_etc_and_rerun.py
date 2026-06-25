#!/usr/bin/env python3
"""Apply ETC fixes to fsp237_minimal_glucose:
  1. Flip stoichiometry of frxn08975_c0 to canonical NDE direction.
  2. Block frxn08975_c0 (lb=ub=0) so it cannot carry flux even if anyone re-opens it.
  3. Set NGAM (non-growth ATP maintenance) on rxn00062_c0 to force ATP demand
     beyond what substrate-level phosphorylation alone can supply. Without this
     the biomass-only objective has degenerate optima and picks 0 flux through
     Complex III/IV. Standard yeast iMM904 NGAM = 1.0 mmol/gDW/hr.
Then re-optimize, report ETC component fluxes, save model+TSV+repainted Escher map.

Run with:  /opt/env/modelseed/bin/python fix_etc_and_rerun.py
"""
import json
import os
import re as _re
import sys

import cobra
from escher import Builder

BASE = '/home/janakae/fungalTemplate/imm904CobraModel'
OUT  = f'{BASE}/fsp237_minimal_glucose'
MODEL_JSON = f'{OUT}/fsp237_minimal_glucose.json'
FLUX_TSV   = f'{OUT}/fsp237_minimal_glucose_fluxes.tsv'
MAP_HTML   = f'{OUT}/fsp237_minimal_glucose_flux_map.html'
ESCHER_MAP = f'{BASE}/iMM904_Central_carbon_metabolism_March28.json'

NGAM_RXN = 'rxn00062_c0'   # ATP + H2O -> ADP + Pi + H+ (cyto)
NGAM_VALUE = 1.0           # standard yeast iMM904 NGAM, mmol/gDW/hr

ETC_RXNS = {
    'rxn09523_m0': 'Complex I (mito NADH dehydrog.)',
    'rxn09524_c0': 'NADH dehydrog. cyto/mito (NDE)',
    'frxn08975_c0': 'External alt. NADH dehydrog. [duplicate]',
    'rxn09417_m0': 'Complex II / succinate dehydrog.',
    'rxn09522_m0': 'Succinate dehydrog. (ubiquinone)',
    'frxn35348_m0': 'Complex III (ubiquinol-cyt c oxidoreductase)',
    'rxn35347_m0': 'Complex IV (cytochrome c oxidase)',
    'rxn08173_m0': 'F(1)-ATPase / ATP synthase',
    'rxn09563_c0': 'Dihydroorotate dehydrog. (ubiquinol producer)',
    NGAM_RXN:      'NGAM (ATP maintenance, cyto)',
}

def report_etc(model, sol, header):
    print(f'\n========== {header} ==========')
    print(f'biomass: {sol.objective_value:.6g}')
    print(f'{"ID":<18} {"flux":>14}   reaction')
    for rid, label in ETC_RXNS.items():
        if rid in [r.id for r in model.reactions]:
            f = sol.fluxes.get(rid, 0)
            mark = '   ' if abs(f) > 1e-9 else ' . '
            print(f'{mark}{rid:<15} {f:>+14.6g}   {label}')
        else:
            print(f' ? {rid:<15} {"(missing)":>14}   {label}')

def main():
    print(f'loading {MODEL_JSON}')
    model = cobra.io.load_json_model(MODEL_JSON)

    # ---- baseline ----
    sol0 = model.optimize()
    report_etc(model, sol0, 'BEFORE fix (current published state)')

    # ---- fix 1+2: flip stoichiometry (if needed), block bounds ----
    # canonical NDE direction: NADH + H+ + ubiquinone -> NAD + ubiquinol
    #   substrates contain cpd00004 (NADH), products contain cpd00003 (NAD)
    r = model.reactions.get_by_id('frxn08975_c0')
    print(f'\nfrxn08975_c0 BEFORE:')
    print(f'  bounds: ({r.lower_bound}, {r.upper_bound})')
    print(f'  reaction: {r.build_reaction_string(use_metabolite_names=False)}')

    nadh_coef = next((c for m, c in r.metabolites.items() if m.id == 'cpd00004_c0'), 0)
    if nadh_coef > 0:
        # NADH is on the product side -> wrong direction; flip
        new_stoich = {met: -coef for met, coef in r.metabolites.items()}
        r.subtract_metabolites(r.metabolites)
        r.add_metabolites(new_stoich)
        print('  -> stoichiometry FLIPPED to canonical NDE direction')
    else:
        print('  -> stoichiometry already canonical; no flip needed')
    r.bounds = (0.0, 0.0)

    print(f'frxn08975_c0 AFTER:')
    print(f'  bounds: ({r.lower_bound}, {r.upper_bound})  [BLOCKED]')
    print(f'  reaction (corrected direction, but disabled): '
          f'{r.build_reaction_string(use_metabolite_names=False)}')

    # ---- fix 3: set NGAM ----
    ngam = model.reactions.get_by_id(NGAM_RXN)
    print(f'\n{NGAM_RXN} BEFORE: bounds ({ngam.lower_bound}, {ngam.upper_bound})')
    ngam.lower_bound = NGAM_VALUE
    print(f'{NGAM_RXN} AFTER:  bounds ({ngam.lower_bound}, {ngam.upper_bound})  '
          f'[NGAM forced ≥ {NGAM_VALUE} mmol/gDW/hr]')

    # ---- re-optimize ----
    sol1 = model.optimize()
    if sol1.objective_value is None:
        print('\nERROR: infeasible after fix; reverting NGAM and aborting.')
        ngam.lower_bound = 0.0
        sys.exit(1)
    report_etc(model, sol1, 'AFTER fix (blocked + corrected + NGAM)')

    if sol1.objective_value is None or sol1.objective_value < 1e-6:
        print('\nERROR: model infeasible after fix; aborting save.')
        sys.exit(1)

    # ---- check ETC functionality ----
    etc_active = {rid: sol1.fluxes.get(rid, 0) for rid in ETC_RXNS
                  if rid in [r.id for r in model.reactions]
                  and abs(sol1.fluxes.get(rid, 0)) > 1e-9}
    complex_iii_active = abs(sol1.fluxes.get('frxn35348_m0', 0)) > 1e-9
    complex_iv_active  = abs(sol1.fluxes.get('rxn35347_m0', 0)) > 1e-9
    nde_active = (abs(sol1.fluxes.get('rxn09523_m0', 0)) > 1e-9
                  or abs(sol1.fluxes.get('rxn09524_c0', 0)) > 1e-9)
    print(f'\nETC summary:')
    print(f'  Complex III (frxn35348_m0) active: {complex_iii_active}')
    print(f'  Complex IV  (rxn35347_m0 ) active: {complex_iv_active}')
    print(f'  Complex I-like NADH dehydrog. active: {nde_active}')

    # ---- save model + fluxes ----
    cobra.io.save_json_model(model, MODEL_JSON)
    sol1.fluxes.to_csv(FLUX_TSV, sep='\t', header=['flux'])
    print(f'\nsaved: {MODEL_JSON}')
    print(f'saved: {FLUX_TSV}')

    # ---- repaint Escher map (same color scheme as build notebook) ----
    flux_dict = sol1.fluxes[sol1.fluxes.abs() > 1e-6].to_dict()
    print(f'active rxns (|flux|>1e-6): {len(flux_dict)}')

    with open(ESCHER_MAP) as f:
        map_json_str = f.read()
    map_data = json.loads(map_json_str)
    map_rxn_ids = {r['bigg_id'] for r in map_data[1]['reactions'].values()}
    print(f'active on central-carbon map: '
          f'{sum(1 for r in flux_dict if r in map_rxn_ids)} / {len(map_rxn_ids)}')

    COMPARTMENT_COLORS = {
        'c0': ('#ff69b4', '#ff1493'),
        'r0': ('#aec7e8', '#1f77b4'),
        'm0': ('#98df8a', '#2ca02c'),
        'x0': ('#ffbb78', '#ff7f0e'),
        'e0': ('#c49c94', '#8c564b'),
        'n0': ('#9edae5', '#17becf'),
        'g0': ('#dbdb8d', '#bcbd22'),
        'v0': ('#c5b0d5', '#9467bd'),
    }
    DEFAULT_COLORS = ('#cc0000', '#cc0000')
    IMM_ESSENTIAL_RED = '#cc0000'
    IMM_ESSENTIALS = {'rxn08475_c0', 'rxn00881_c0', 'ORNt3m_c0'}

    def _hex_interp(low, high, t):
        h2r = lambda h: tuple(int(h.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
        r1, g1, b1 = h2r(low); r2, g2, b2 = h2r(high)
        return '#{:02x}{:02x}{:02x}'.format(
            int(r1 + (r2 - r1) * t),
            int(g1 + (g2 - g1) * t),
            int(b1 + (b2 - b1) * t),
        )

    max_abs = max((abs(v) for v in flux_dict.values()), default=1.0) or 1.0
    color_map, size_map = {}, {}
    for rid, flux in flux_dict.items():
        if abs(flux) <= 1e-9: continue
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue
        suffixes = [met.id.rsplit('_', 1)[-1] for met in rxn.metabolites if '_' in met.id]
        comp = max(set(suffixes), key=suffixes.count) if suffixes else 'c0'
        t = 0.35 + 0.65 * (abs(flux) / max_abs) ** 0.35
        if rid in IMM_ESSENTIALS:
            color_map[rid] = IMM_ESSENTIAL_RED
        else:
            low, high = COMPARTMENT_COLORS.get(comp, DEFAULT_COLORS)
            color_map[rid] = _hex_interp(low, high, t)
        size_map[rid] = 3 + t * 18

    builder = Builder(
        map_json=map_json_str,
        reaction_data=flux_dict,
        reaction_scale=[
            {'type': 'min',    'color': '#cccccc', 'size': 2},
            {'type': 'value',  'value': 0, 'color': '#cccccc', 'size': 2},
            {'type': 'median', 'color': '#888888', 'size': 10},
            {'type': 'max',    'color': '#111111', 'size': 20},
        ],
        reaction_no_data_color='#dcdcdc',
        reaction_no_data_size=3,
    )
    builder.save_html(MAP_HTML)

    with open(MAP_HTML) as f:
        html = f.read()
    html = _re.sub(r'<script>\s*\(function\(\).*?\}\)\(\);\s*</script>', '', html, flags=_re.DOTALL)
    html = _re.sub(r'<div id="compartment-legend".*?</div>\s*', '', html, flags=_re.DOTALL)

    color_js = json.dumps(color_map); size_js = json.dumps(size_map)
    legend_items = [
        ('#ff1493', 'Cytosol (c0)'),
        ('#1f77b4', 'ER (r0)'),
        ('#2ca02c', 'Mitochondria (m0)'),
        ('#ff7f0e', 'Peroxisome (x0)'),
        ('#8c564b', 'Extracellular (e0)'),
        ('#17becf', 'Nucleus (n0)'),
        ('#bcbd22', 'Golgi (g0)'),
        ('#9467bd', 'Vacuole (v0)'),
        ('#cc0000', 'iMM904 essentials (kept)'),
        ('#dcdcdc', 'No flux / no data'),
    ]
    legend_html = ''.join(
        f'<div style="display:flex;align-items:center;margin-bottom:5px">'
        f'<div style="width:26px;height:4px;background:{c};margin-right:8px;border-radius:2px"></div>'
        f'<span>{lbl}</span></div>'
        for c, lbl in legend_items
    )
    g = sol1.objective_value
    inject = f"""
<script>
(function() {{
  var colorMap = {color_js};
  var sizeMap  = {size_js};
  function applyColors() {{
    var applied = 0;
    document.querySelectorAll('g.reaction').forEach(function(g) {{
      var labelEl = g.querySelector('.reaction-label');
      if (!labelEl) return;
      var rxnId = labelEl.textContent.trim().split(/\\s+/)[0];
      if (colorMap[rxnId]) {{
        var color = colorMap[rxnId];
        var size  = (sizeMap[rxnId] || 2) + 'px';
        g.querySelectorAll('path, line, polyline').forEach(function(p) {{
          p.style.stroke = color;
          p.style.strokeWidth = size;
        }});
        labelEl.style.fill = color;
        applied++;
      }}
    }});
    console.log('Compartment colors applied to', applied, 'active reactions');
  }}
  var attempts = 0;
  var interval = setInterval(function() {{
    if (document.querySelectorAll('g.reaction').length > 0) {{
      clearInterval(interval);
      applyColors();
    }} else if (++attempts > 60) {{
      clearInterval(interval);
    }}
  }}, 200);
}})();
</script>

<div id="compartment-legend" style="
  position:fixed; bottom:20px; right:20px;
  background:rgba(255,255,255,0.97);
  border:1px solid #ddd; border-radius:10px;
  padding:16px 20px; font-family:sans-serif;
  font-size:12px; z-index:9999;
  box-shadow:0 3px 12px rgba(0,0,0,0.18)">
  <b style="font-size:13px">fsp237 minimal model &mdash; ETC fix applied</b>
  <div style="margin:8px 0">{legend_html}</div>
  <i style="color:#999;font-size:11px">Intensity &amp; thickness = flux magnitude</i><br>
  <i style="color:#999;font-size:11px">frxn08975_c0 blocked &middot; biomass={g:.4f}</i>
</div>
"""
    html = html.replace('</body>', inject + '\n</body>')
    with open(MAP_HTML, 'w') as f:
        f.write(html)
    print(f'wrote {MAP_HTML} ({os.path.getsize(MAP_HTML)/1024:.1f} KB)')

if __name__ == '__main__':
    main()
