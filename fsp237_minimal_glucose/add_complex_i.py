#!/usr/bin/env python3
"""Add proton-pumping Complex I (frxn13726_m0) to fsp237_minimal_glucose.

Reaction (user-supplied):
  NADH[m0] + Ubiquinone-6[m0] + 5 H+[m0]  =>  NAD[m0] + Ubiquinol-6[m0] + 4 H+[c0]

Idempotent: if frxn13726_m0 already exists, only verifies bounds.
After adding, re-optimizes and reports the full ETC, then saves model + TSV
+ repainted Escher map.
"""
import json
import os
import re as _re
import sys

import cobra
from cobra import Reaction
from escher import Builder

BASE = '/home/janakae/fungalTemplate/imm904CobraModel'
OUT  = f'{BASE}/fsp237_minimal_glucose'
MODEL_JSON = f'{OUT}/fsp237_minimal_glucose.json'
FLUX_TSV   = f'{OUT}/fsp237_minimal_glucose_fluxes.tsv'
MAP_HTML   = f'{OUT}/fsp237_minimal_glucose_flux_map.html'
ESCHER_MAP = f'{BASE}/iMM904_Central_carbon_metabolism_March28.json'

NEW_RXN_ID = 'frxn13726_m0'
NEW_RXN_NAME = 'NADH:ubiquinone oxidoreductase, proton-pumping (Complex I)_m0'
# substrates negative, products positive
NEW_RXN_STOICH = {
    'cpd00004_m0': -1,   # NADH
    'cpd15290_m0': -1,   # Ubiquinone-6
    'cpd00067_m0': -5,   # 5 H+ matrix consumed (1 chemistry + 4 pumped)
    'cpd00003_m0':  1,   # NAD
    'cpd15291_m0':  1,   # Ubiquinol-6
    'cpd00067_c0':  4,   # 4 H+ pumped to IMS (modeled as c0)
}
NEW_RXN_BOUNDS = (0.0, 1000.0)   # forward only

ETC_RXNS = {
    NEW_RXN_ID:    'Complex I, proton-pumping (NEW)',
    'rxn09523_m0': 'Complex I dehydrogenase (NDI1-like)',
    'rxn09524_c0': 'NADH dehydrog. cyto/mito (NDE-like)',
    'frxn08975_c0': 'External alt. NADH dehydrog. [blocked]',
    'rxn09417_m0': 'Complex II / succinate dehydrog.',
    'rxn09522_m0': 'Succinate dehydrog. (ubiquinone)',
    'frxn35348_m0': 'Complex III (ubiquinol-cyt c oxidoreductase)',
    'rxn35347_m0': 'Complex IV (cytochrome c oxidase)',
    'rxn08173_m0': 'F(1)-ATPase / ATP synthase',
    'rxn00062_c0': 'NGAM (ATP maintenance, cyto)',
    'rxn09563_c0': 'Dihydroorotate dehydrog. (ubiquinol producer)',
    'EX_cpd00007_e0': 'O2 exchange',
    'EX_cpd00027_e0': 'Glucose exchange',
    'bio1': 'Biomass',
}

def report(model, sol, header):
    print(f'\n========== {header} ==========')
    print(f'biomass: {sol.objective_value:.6g}')
    rxn_ids = {r.id for r in model.reactions}
    print(f'{"ID":<20} {"flux":>14}   reaction')
    for rid, label in ETC_RXNS.items():
        if rid in rxn_ids:
            f = sol.fluxes.get(rid, 0)
            mark = '   ' if abs(f) > 1e-9 else ' . '
            print(f'{mark}{rid:<17} {f:>+14.6g}   {label}')
        else:
            print(f' ? {rid:<17} {"(missing)":>14}   {label}')

def main():
    print(f'loading {MODEL_JSON}')
    model = cobra.io.load_json_model(MODEL_JSON)

    sol0 = model.optimize()
    report(model, sol0, 'BEFORE adding Complex I (proton-pumping)')

    # ---- add or update frxn13726_m0 ----
    if NEW_RXN_ID in [r.id for r in model.reactions]:
        print(f'\n{NEW_RXN_ID} already exists; updating stoichiometry + bounds')
        r = model.reactions.get_by_id(NEW_RXN_ID)
        r.subtract_metabolites(r.metabolites)
    else:
        print(f'\nadding new reaction {NEW_RXN_ID}')
        r = Reaction(NEW_RXN_ID)
        r.name = NEW_RXN_NAME
        model.add_reactions([r])

    # Build the stoich using the model's existing metabolite objects
    met_dict = {}
    for mid, coef in NEW_RXN_STOICH.items():
        try:
            met_dict[model.metabolites.get_by_id(mid)] = coef
        except KeyError:
            sys.exit(f'ERROR: metabolite {mid} not in model')
    r.add_metabolites(met_dict)
    r.bounds = NEW_RXN_BOUNDS
    r.name = NEW_RXN_NAME
    print(f'  bounds: {r.bounds}')
    print(f'  eq: {r.build_reaction_string(use_metabolite_names=False)}')

    # ---- re-optimize ----
    sol1 = model.optimize()
    if sol1.objective_value is None:
        sys.exit('\nERROR: infeasible after adding Complex I; aborting save.')
    report(model, sol1, 'AFTER adding Complex I (proton-pumping)')

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
  <b style="font-size:13px">fsp237 + Complex I (proton-pumping) added</b>
  <div style="margin:8px 0">{legend_html}</div>
  <i style="color:#999;font-size:11px">Intensity &amp; thickness = flux magnitude</i><br>
  <i style="color:#999;font-size:11px">frxn13726_m0 added &middot; biomass={g:.4f}</i>
</div>
"""
    html = html.replace('</body>', inject + '\n</body>')
    with open(MAP_HTML, 'w') as f:
        f.write(html)
    print(f'wrote {MAP_HTML} ({os.path.getsize(MAP_HTML)/1024:.1f} KB)')

if __name__ == '__main__':
    main()
