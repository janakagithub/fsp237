#!/opt/env/modelseed/bin/python3
"""Build an Escher map coloured by S1 expression per reaction (analogous to
the flux-painted map but with |flux| replaced by aggregated expression score).

Re-uses the exact same compartment palette + intensity ramp + legend layout
as build_escher_maps.py, so the expression map sits alongside the aerobic /
anaerobic flux maps as a third, directly comparable overlay on the same
iMM904 central-carbon layout.

  intensity = 0.15 + 0.85 * (expr / P95) ^ 0.60
              (P95 = 95th-percentile score across scored reactions)

Reactions with no expression score (no GPR, or all GPR genes missing) draw
in the standard 'no data' grey.

Output: /home/janakae/fsp237/atp-safe/map_expression.html
"""
import json
import os
import re
import sys
from pathlib import Path

import cobra
import escher
import pandas as pd

BASE = Path('/home/janakae/fungalTemplate/imm904CobraModel')
MODEL_PATH = BASE / 'simulations/gapfill_v1_v2/models/fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json'
MAP_JSON = BASE / 'iMM904_Central_carbon_metabolism_March28.json'
REXP_TSV = BASE / 'rna_seq_integration/outputs/reaction_expression.tsv'
OUT_HTML = Path('/home/janakae/fsp237/atp-safe/map_expression.html')

# Same palette as build_escher_maps.py
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


def hex_interp(low, high, t):
    h2r = lambda h: tuple(int(h.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
    r1, g1, b1 = h2r(low); r2, g2, b2 = h2r(high)
    return '#{:02x}{:02x}{:02x}'.format(
        int(r1 + (r2 - r1) * t),
        int(g1 + (g2 - g1) * t),
        int(b1 + (b2 - b1) * t),
    )


def compartment_for(rxn):
    sufs = [met.id.rsplit('_', 1)[-1] for met in rxn.metabolites if '_' in met.id]
    return max(set(sufs), key=sufs.count) if sufs else 'c0'


def main():
    print(f'model      : {MODEL_PATH.name}')
    print(f'expression : {REXP_TSV.name}')
    print(f'map        : {MAP_JSON.name}')
    model = cobra.io.load_json_model(str(MODEL_PATH))
    rexp = pd.read_csv(REXP_TSV, sep='\t', comment='#').set_index('rxn_id')
    expr_by_rxn = rexp['agg_mean_log2TPMp1'].to_dict()

    # Numeric score dict for reactions that have a numeric expression
    scored = {rid: float(v) for rid, v in expr_by_rxn.items()
              if v is not None and pd.notna(v)}
    if not scored:
        raise SystemExit('no reactions with numeric expression scores')
    vals = sorted(scored.values())
    p95 = vals[int(0.95 * len(vals))]
    print(f'scored reactions: {len(scored)}   P95 score: {p95:.3f}   max: {max(vals):.3f}')

    color_map, size_map = {}, {}
    for rid, s in scored.items():
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue
        comp = compartment_for(rxn)
        norm = min(1.0, max(0.0, s / p95))
        t = 0.15 + 0.85 * (norm ** 0.60)
        low, high = COMPARTMENT_COLORS.get(comp, DEFAULT_COLORS)
        color_map[rid] = hex_interp(low, high, t)
        size_map[rid] = 3 + t * 18

    with open(MAP_JSON) as f:
        map_json_str = f.read()

    # Escher's own paint is a plain grayscale; JS overlay repaints.
    builder = escher.Builder(
        map_json=map_json_str,
        reaction_data=scored,
        reaction_scale=[
            {'type': 'min',    'color': '#cccccc', 'size': 2},
            {'type': 'value',  'value': 0, 'color': '#cccccc', 'size': 2},
            {'type': 'median', 'color': '#888888', 'size': 10},
            {'type': 'max',    'color': '#111111', 'size': 20},
        ],
        reaction_no_data_color='#e7ebef',
        reaction_no_data_size=3,
    )
    builder.save_html(str(OUT_HTML))

    with open(OUT_HTML) as f:
        html = f.read()
    html = re.sub(r'<script id="atp-safe-recolor">.*?</script>', '', html, flags=re.DOTALL)
    html = re.sub(r'<div id="atp-safe-legend".*?</div>\s*(?=</body>)', '', html, flags=re.DOTALL)

    color_js = json.dumps(color_map)
    size_js = json.dumps(size_map)

    legend_items = [
        ('#ff1493', 'Cytosol (c0)'),
        ('#1f77b4', 'ER (r0)'),
        ('#2ca02c', 'Mitochondria (m0)'),
        ('#ff7f0e', 'Peroxisome (x0)'),
        ('#8c564b', 'Extracellular (e0)'),
        ('#17becf', 'Nucleus (n0)'),
        ('#bcbd22', 'Golgi (g0)'),
        ('#9467bd', 'Vacuole (v0)'),
        ('#e7ebef', 'No expression / no GPR'),
    ]
    legend_html = ''.join(
        f'<div style="display:flex;align-items:center;margin-bottom:5px">'
        f'<div style="width:26px;height:4px;background:{c};margin-right:8px;border-radius:2px"></div>'
        f'<span>{lbl}</span></div>'
        for c, lbl in legend_items
    )

    n_scored = len(scored)
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
    console.log('atp-safe expression recolor: painted', applied, 'reactions');
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
  <b style="font-size:13px">fsp237 ATP-safe GSM &mdash; S1 expression</b>
  <div style="margin:8px 0">{legend_html}</div>
  <i style="color:#999;font-size:11px">Intensity &amp; thickness = mean log<sub>2</sub>(TPM+1) aggregated over GPR</i><br>
  <i style="color:#999;font-size:11px">{n_scored} scored reactions &middot; P95 = {p95:.2f}</i>
</div>
"""

    html = html.replace('</body>', inject + '\n</body>')
    with open(OUT_HTML, 'w') as f:
        f.write(html)
    print(f'wrote {OUT_HTML}  ({os.path.getsize(OUT_HTML)/1024:.0f} KB, {len(color_map)} reactions recoloured)')


if __name__ == '__main__':
    main()
