#!/opt/env/modelseed/bin/python3
"""Build an Escher map painted by proteomics-derived GIMME flux (PDA baseline),
the proteome analog of map_expression.html / map_aerobic.html. Same iMM904
central-carbon layout, same compartment palette + intensity ramp + JS recolor,
so it sits directly alongside the transcriptome and flux maps.

intensity = 0.15 + 0.85 * (|flux| / P95) ^ 0.60

Output: atp-safe/map_proteomics.html
"""
import json
import os
import re
from pathlib import Path

import cobra
import escher
import pandas as pd

ROOT = Path("/home/janakae/fsp237")
MODEL_PATH = (ROOT / "simulations/gapfill_v1_v2/models/"
              "fsp237_gapfilled_Version10_vlcfa_complete_genes_integrated.json")
MAP_JSON = Path("/home/janakae/fungalTemplate/imm904CobraModel/"
                "iMM904_Central_carbon_metabolism_March28.json")
GIMME_FM = ROOT / "proteomics_integration/outputs/context/gimme_flux_matrix.tsv"
OUT_HTML = ROOT / "atp-safe/map_proteomics.html"

COMPARTMENT_COLORS = {
    'c0': ('#ff69b4', '#ff1493'), 'r0': ('#aec7e8', '#1f77b4'),
    'm0': ('#98df8a', '#2ca02c'), 'x0': ('#ffbb78', '#ff7f0e'),
    'e0': ('#c49c94', '#8c564b'), 'n0': ('#9edae5', '#17becf'),
    'g0': ('#dbdb8d', '#bcbd22'), 'v0': ('#c5b0d5', '#9467bd'),
}
DEFAULT_COLORS = ('#cc0000', '#cc0000')


def hex_interp(low, high, t):
    h2r = lambda h: tuple(int(h.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
    r1, g1, b1 = h2r(low); r2, g2, b2 = h2r(high)
    return '#{:02x}{:02x}{:02x}'.format(int(r1 + (r2 - r1) * t),
                                        int(g1 + (g2 - g1) * t),
                                        int(b1 + (b2 - b1) * t))


def compartment_for(rxn):
    sufs = [m.id.rsplit('_', 1)[-1] for m in rxn.metabolites if '_' in m.id]
    return max(set(sufs), key=sufs.count) if sufs else 'c0'


def main():
    model = cobra.io.load_json_model(str(MODEL_PATH))
    gf = pd.read_csv(GIMME_FM, sep="\t", index_col=0)
    flux = {rid: abs(float(v)) for rid, v in gf["flux_PDA"].dropna().items()
            if abs(float(v)) > 1e-6}
    if not flux:
        raise SystemExit("no GIMME PDA flux to paint")
    vals = sorted(flux.values())
    p95 = vals[int(0.95 * len(vals))]
    print(f"painted reactions: {len(flux)}  P95 |flux|: {p95:.3f}  max: {max(vals):.3f}")

    color_map, size_map = {}, {}
    for rid, s in flux.items():
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue
        norm = min(1.0, max(0.0, s / p95))
        t = 0.15 + 0.85 * (norm ** 0.60)
        low, high = COMPARTMENT_COLORS.get(compartment_for(rxn), DEFAULT_COLORS)
        color_map[rid] = hex_interp(low, high, t)
        size_map[rid] = 3 + t * 18

    builder = escher.Builder(
        map_json=open(MAP_JSON).read(), reaction_data=flux,
        reaction_scale=[{'type': 'min', 'color': '#cccccc', 'size': 2},
                        {'type': 'value', 'value': 0, 'color': '#cccccc', 'size': 2},
                        {'type': 'median', 'color': '#888888', 'size': 10},
                        {'type': 'max', 'color': '#111111', 'size': 20}],
        reaction_no_data_color='#e7ebef', reaction_no_data_size=3)
    builder.save_html(str(OUT_HTML))

    html = open(OUT_HTML).read()
    html = re.sub(r'<script id="atp-safe-recolor">.*?</script>', '', html, flags=re.DOTALL)
    html = re.sub(r'<div id="atp-safe-legend".*?</div>\s*(?=</body>)', '', html, flags=re.DOTALL)

    legend_items = [('#ff1493', 'Cytosol (c0)'), ('#1f77b4', 'ER (r0)'),
                    ('#2ca02c', 'Mitochondria (m0)'), ('#ff7f0e', 'Peroxisome (x0)'),
                    ('#8c564b', 'Extracellular (e0)'), ('#17becf', 'Nucleus (n0)'),
                    ('#bcbd22', 'Golgi (g0)'), ('#9467bd', 'Vacuole (v0)'),
                    ('#e7ebef', 'No flux / no GPR')]
    legend_html = ''.join(
        f'<div style="display:flex;align-items:center;margin-bottom:5px">'
        f'<div style="width:26px;height:4px;background:{c};margin-right:8px;border-radius:2px"></div>'
        f'<span>{lbl}</span></div>' for c, lbl in legend_items)

    inject = f"""
<script id="atp-safe-recolor">
(function() {{
  var colorMap = {json.dumps(color_map)};
  var sizeMap  = {json.dumps(size_map)};
  function applyColors() {{
    var applied = 0;
    document.querySelectorAll('g.reaction').forEach(function(g) {{
      var labelEl = g.querySelector('.reaction-label');
      if (!labelEl) return;
      var rxnId = labelEl.textContent.trim().split(/\\s+/)[0];
      if (colorMap[rxnId]) {{
        var color = colorMap[rxnId], size = (sizeMap[rxnId] || 2) + 'px';
        g.querySelectorAll('path, line, polyline').forEach(function(p) {{
          p.style.stroke = color; p.style.strokeWidth = size; }});
        labelEl.style.fill = color; applied++;
      }}
    }});
    console.log('atp-safe proteomics recolor: painted', applied, 'reactions');
  }}
  var attempts = 0;
  var interval = setInterval(function() {{
    if (document.querySelectorAll('g.reaction').length > 0) {{
      clearInterval(interval); applyColors();
    }} else if (++attempts > 60) {{ clearInterval(interval); }}
  }}, 200);
}})();
</script>
<div id="atp-safe-legend" style="position:fixed; bottom:20px; right:20px;
  background:rgba(255,255,255,0.97); border:1px solid #ddd; border-radius:10px;
  padding:16px 20px; font-family:sans-serif; font-size:12px; z-index:9999;
  box-shadow:0 3px 12px rgba(0,0,0,0.18); max-width:280px">
  <b style="font-size:13px">fsp237 &mdash; proteome-constrained flux (GIMME, full PDB)</b>
  <div style="margin:8px 0">{legend_html}</div>
  <i style="color:#999;font-size:11px">Intensity &amp; thickness = |GIMME flux| from the measured proteome</i><br>
  <i style="color:#999;font-size:11px">{len(color_map)} reactions painted &middot; P95 = {p95:.2f}</i>
</div>
"""
    html = html.replace('</body>', inject + '\n</body>')
    open(OUT_HTML, "w").write(html)
    print(f"wrote {OUT_HTML} ({os.path.getsize(OUT_HTML)/1024:.0f} KB, "
          f"{len(color_map)} reactions painted)")


if __name__ == "__main__":
    main()
