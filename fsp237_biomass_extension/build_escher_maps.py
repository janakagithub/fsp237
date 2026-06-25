#!/usr/bin/env python3
"""Regenerate Escher map HTMLs for /home/janakae/fsp237/atp-safe/ with
fluxes from fsp237_atp_safe_gsm_extended.json, styled the same way as
fsp237_minimal_glucose/BuildMinimalFSP237Model.ipynb (cell 19):

  - call Escher Builder with a plain grayscale `reaction_scale` so Escher
    paints the map with the standard min/median/max palette
  - then post-process the saved HTML, injecting a small browser-side script
    that recolors each <g.reaction> using:
        * the reaction's primary compartment  -> base palette (cytosol pink,
          mito green, ER blue, etc.)
        * intensity = 0.35 + 0.65 * (|flux|/max_abs)^0.35  -> hex-interpolated
          between low/high stops of that compartment
        * stroke width scaled the same way
        * (optional) any reaction in HIGHLIGHT_SET gets a bright red override
    plus a floating bottom-right legend with palette + biomass info.

This is the same look used on the minimal-glucose flux map at
fsp237/index.html ('original GSM with ETC fixes'); re-applying it here so the
atp-safe maps match that visual style consistently across the site.
"""
import json
import os
import re

import cobra
import escher

BASE = '/home/janakae/fungalTemplate/imm904CobraModel'
# Source the LATEST V6 model (gapfilled + dedup + dirlock + gene-integrated);
# see SUMMARY.md in simulations/gapfill_v1_v2/reports/ for the version trail.
MODEL_PATH = f'{BASE}/simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json'
GSM_PATH   = f'{BASE}/fsp237_minimal_glucose/fsp237_minimal_glucose.json'
MAP_JSON   = f'{BASE}/iMM904_Central_carbon_metabolism_March28.json'

OUT_DIR = '/home/janakae/fsp237/atp-safe'
OUT_AER = f'{OUT_DIR}/map_aerobic.html'
OUT_ANA = f'{OUT_DIR}/map_anaerobic.html'

# Palette: (low, high) per compartment -- same as the minimal-glucose notebook.
COMPARTMENT_COLORS = {
    'c0': ('#ff69b4', '#ff1493'),   # cytosol: pink
    'r0': ('#aec7e8', '#1f77b4'),   # ER: blue
    'm0': ('#98df8a', '#2ca02c'),   # mitochondria: green
    'x0': ('#ffbb78', '#ff7f0e'),   # peroxisome: orange
    'e0': ('#c49c94', '#8c564b'),   # extracellular: brown
    'n0': ('#9edae5', '#17becf'),   # nucleus: teal
    'g0': ('#dbdb8d', '#bcbd22'),   # Golgi: olive
    'v0': ('#c5b0d5', '#9467bd'),   # vacuole: purple
}
DEFAULT_COLORS  = ('#cc0000', '#cc0000')
HIGHLIGHT_COLOR = '#cc0000'

# Reactions added in the biomass-extension pass (cell-wall, mannitol, melanin,
# arabinose, the new alpha-1,3-glucan synthase). Highlighted bright red on
# the map when they happen to lie on the iMM904 CCM layout. Read from the
# extension report so this list stays in sync with extend_biomass.py.
EXTENSION_REPORT = f'{BASE}/fsp237_biomass_extension/extension_report.tsv'


def load_highlight_set():
    if not os.path.exists(EXTENSION_REPORT):
        return set()
    out = set()
    with open(EXTENSION_REPORT) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        try:
            i_id, i_st = header.index('rxn_id'), header.index('status')
        except ValueError:
            return set()
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) > max(i_id, i_st) and parts[i_st] == 'added':
                out.add(parts[i_id])
    return out


def hex_interp(low, high, t):
    h2r = lambda h: tuple(int(h.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
    r1, g1, b1 = h2r(low); r2, g2, b2 = h2r(high)
    return '#{:02x}{:02x}{:02x}'.format(
        int(r1 + (r2 - r1) * t),
        int(g1 + (g2 - g1) * t),
        int(b1 + (b2 - b1) * t),
    )


def fba_with_gsm_media(model, gsm, objective_id, aerobic=True):
    """Set GSM-style media + O2, optimize, return (obj_value, signed_flux_dict)."""
    with model:
        gsm_bounds = {r.id: (r.lower_bound, r.upper_bound) for r in gsm.exchanges}
        for ex in model.exchanges:
            if ex.id in gsm_bounds:
                ex.lower_bound, ex.upper_bound = gsm_bounds[ex.id]
            else:
                ex.lower_bound = 0
        if not aerobic and 'EX_cpd00007_e0' in [r.id for r in model.reactions]:
            model.reactions.get_by_id('EX_cpd00007_e0').lower_bound = 0
        for r in model.reactions:
            r.objective_coefficient = 1 if r.id == objective_id else 0
        sol = model.optimize()
        return (sol.objective_value or 0.0), {rid: float(v) for rid, v in sol.fluxes.items()}


def compartment_for(rxn):
    sufs = [met.id.rsplit('_', 1)[-1] for met in rxn.metabolites if '_' in met.id]
    return max(set(sufs), key=sufs.count) if sufs else 'c0'


def build_map(model, map_json_str, flux_dict, out_path, biomass_val, atp_yield,
              condition_label, highlight_set):
    """Build one Escher map HTML with compartment-coloured per-reaction overlays."""
    # Trim near-zero noise (Escher would draw them anyway as min color)
    active = {k: v for k, v in flux_dict.items() if abs(v) > 1e-9}
    max_abs = max((abs(v) for v in active.values()), default=1.0) or 1.0

    color_map, size_map = {}, {}
    for rid, flux in active.items():
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue
        comp = compartment_for(rxn)
        t = 0.35 + 0.65 * (abs(flux) / max_abs) ** 0.35
        if rid in highlight_set:
            color_map[rid] = HIGHLIGHT_COLOR
        else:
            low, high = COMPARTMENT_COLORS.get(comp, DEFAULT_COLORS)
            color_map[rid] = hex_interp(low, high, t)
        size_map[rid] = 3 + t * 18

    # Base Escher render -- simple grayscale ramp; the post-process script
    # paints over it. (Keeping Escher's own paint in place means the map is
    # still useful even if the inline script hasn't run yet.)
    builder = escher.Builder(
        map_json=map_json_str,
        reaction_data=active,
        reaction_scale=[
            {'type': 'min',    'color': '#cccccc', 'size': 2},
            {'type': 'value',  'value': 0, 'color': '#cccccc', 'size': 2},
            {'type': 'median', 'color': '#888888', 'size': 10},
            {'type': 'max',    'color': '#111111', 'size': 20},
        ],
        reaction_no_data_color='#dcdcdc',
        reaction_no_data_size=3,
    )
    builder.save_html(out_path)

    # Post-process: strip any previously-injected block, then add the script
    # + floating legend.
    with open(out_path) as f:
        html = f.read()
    # Remove any prior injected re-color block (idempotent rebuilds)
    html = re.sub(r'<script id="atp-safe-recolor">.*?</script>', '', html, flags=re.DOTALL)
    html = re.sub(r'<div id="atp-safe-legend".*?</div>\s*(?=</body>)', '', html, flags=re.DOTALL)

    color_js = json.dumps(color_map)
    size_js  = json.dumps(size_map)

    legend_items = [
        ('#ff1493', 'Cytosol (c0)'),
        ('#1f77b4', 'ER (r0)'),
        ('#2ca02c', 'Mitochondria (m0)'),
        ('#ff7f0e', 'Peroxisome (x0)'),
        ('#8c564b', 'Extracellular (e0)'),
        ('#17becf', 'Nucleus (n0)'),
        ('#bcbd22', 'Golgi (g0)'),
        ('#9467bd', 'Vacuole (v0)'),
        ('#cc0000', 'Biomass-extension additions'),
        ('#dcdcdc', 'No flux / no data'),
    ]
    legend_html = ''.join(
        f'<div style="display:flex;align-items:center;margin-bottom:5px">'
        f'<div style="width:26px;height:4px;background:{c};margin-right:8px;border-radius:2px"></div>'
        f'<span>{lbl}</span></div>'
        for c, lbl in legend_items
    )

    inject = f"""
<script id="atp-safe-recolor">
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
    console.log('atp-safe recolor: painted', applied, 'reactions');
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

<div id="atp-safe-legend" style="
  position:fixed; bottom:20px; right:20px;
  background:rgba(255,255,255,0.97);
  border:1px solid #ddd; border-radius:10px;
  padding:16px 20px; font-family:sans-serif;
  font-size:12px; z-index:9999;
  box-shadow:0 3px 12px rgba(0,0,0,0.18); max-width:280px">
  <b style="font-size:13px">fsp237 ATP-safe GSM &mdash; {condition_label}</b>
  <div style="margin:8px 0">{legend_html}</div>
  <i style="color:#999;font-size:11px">Intensity &amp; thickness = flux magnitude</i><br>
  <i style="color:#999;font-size:11px">biomass = {biomass_val:.4f} 1/h &middot; ATP/glc = {atp_yield:.2f}</i>
</div>
"""

    html = html.replace('</body>', inject + '\n</body>')
    with open(out_path, 'w') as f:
        f.write(html)

    n_painted = len(color_map)
    print(f'  wrote {out_path}  ({os.path.getsize(out_path)/1024:.0f} KB; {n_painted} reactions recolored; {condition_label})')


def atp_yield(model, atp_rxn='bio1', aerobic=True, glc=1.0):
    with model:
        for ex in model.exchanges:
            mid = next(iter(ex.metabolites.keys())).id
            if mid in {'cpd00001_e0','cpd00009_e0','cpd00011_e0','cpd00067_e0',
                       'cpd00013_e0','cpd00048_e0','cpd10515_e0','cpd00971_e0',
                       'cpd00205_e0'}:
                ex.lower_bound = -1000
            else:
                ex.lower_bound = 0
        if 'EX_cpd00027_e0' in [r.id for r in model.reactions]:
            model.reactions.get_by_id('EX_cpd00027_e0').lower_bound = -glc
            model.reactions.get_by_id('EX_cpd00027_e0').upper_bound = 0
        if 'EX_cpd00007_e0' in [r.id for r in model.reactions]:
            model.reactions.get_by_id('EX_cpd00007_e0').lower_bound = -1000 if aerobic else 0
        for r in model.reactions:
            r.objective_coefficient = 1 if r.id == atp_rxn else 0
        return model.optimize().objective_value or 0.0


def main():
    print(f'loading model: {MODEL_PATH}')
    model = cobra.io.load_json_model(MODEL_PATH)
    print(f'loading media reference: {GSM_PATH}')
    gsm = cobra.io.load_json_model(GSM_PATH)
    with open(MAP_JSON) as f:
        map_json_str = f.read()
    map_data = json.loads(map_json_str)
    map_rxn_ids = {r['bigg_id'] for r in map_data[1]['reactions'].values()}
    print(f'CCM map reactions: {len(map_rxn_ids)}')

    highlights = load_highlight_set()
    print(f'highlight set (biomass-extension additions): {len(highlights)}')

    print('\n=== aerobic FBA on bio_gsm ===')
    bio_a, flx_a = fba_with_gsm_media(model, gsm, 'bio_gsm', aerobic=True)
    atp_a = atp_yield(model, 'bio1', aerobic=True, glc=1.0)
    print(f'  biomass = {bio_a:.4f} 1/h, ATP/glc = {atp_a:.2f}')
    on_map_a = sum(1 for r in flx_a if r in map_rxn_ids and abs(flx_a[r]) > 1e-9)
    print(f'  active on CCM map: {on_map_a} / {len(map_rxn_ids)}')

    print('\n=== anaerobic FBA on bio_gsm ===')
    bio_n, flx_n = fba_with_gsm_media(model, gsm, 'bio_gsm', aerobic=False)
    atp_n = atp_yield(model, 'bio1', aerobic=False, glc=1.0)
    print(f'  biomass = {bio_n:.4f} 1/h, ATP/glc = {atp_n:.2f}')
    on_map_n = sum(1 for r in flx_n if r in map_rxn_ids and abs(flx_n[r]) > 1e-9)
    print(f'  active on CCM map: {on_map_n} / {len(map_rxn_ids)}')

    print('\n=== building Escher HTMLs ===')
    build_map(model, map_json_str, flx_a, OUT_AER,
              biomass_val=bio_a, atp_yield=atp_a,
              condition_label='aerobic', highlight_set=highlights)
    build_map(model, map_json_str, flx_n, OUT_ANA,
              biomass_val=bio_n, atp_yield=atp_n,
              condition_label='anaerobic', highlight_set=highlights)


if __name__ == '__main__':
    main()
