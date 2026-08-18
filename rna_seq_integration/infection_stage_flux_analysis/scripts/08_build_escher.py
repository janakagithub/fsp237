"""08 - Dynamic Escher flux map for the Infection-Stage Flux tab.

Produces two artifacts under atp-safe/:
  * infection_flux_by_condition.json - compact {method: {condition: {rxn: flux}}}
    (nonzero fluxes only) plus condition labels/stages, built from the wide
    per-method aerobic flux pivots.
  * infection_escher.html - a ONE-OFF Escher builder page cloned from the frozen
    map_aerobic.html (same embedded map_data / model_data), but with the static
    single-condition recolor script and legend REPLACED by a data-driven recolor
    that reads ?method=&cond= from its own URL and repaints from the flux JSON.

index.html embeds this page in an <iframe> and swaps its ?query when the user
picks a method/medium (see loadInfEscher). The frozen map_aerobic/anaerobic/
expression maps are NOT modified.
"""
import os
import re
import json
import pandas as pd
from _common import FLUX_RESULTS, ANALYSIS, METHODS

ATP_SAFE = "/home/janakae/fsp237/atp-safe"
SRC_MAP = os.path.join(ATP_SAFE, "map_aerobic.html")
OUT_MAP = os.path.join(ATP_SAFE, "infection_escher.html")
OUT_JSON = os.path.join(ATP_SAFE, "infection_flux_by_condition.json")
FLUX_EPS = 1e-6

DYNAMIC_RECOLOR = r"""
<script id="infection-recolor">
(function () {
  function esc(s){var d=document.createElement('div');d.textContent=(s==null?'':s);return d.innerHTML;}
  function interp(a,b,t){
    function h(x){return [parseInt(x.slice(1,3),16),parseInt(x.slice(3,5),16),parseInt(x.slice(5,7),16)];}
    var A=h(a),B=h(b),o=[0,1,2].map(function(i){return Math.round(A[i]+(B[i]-A[i])*t);});
    return '#'+o.map(function(v){return v.toString(16).padStart(2,'0');}).join('');
  }
  var q=new URLSearchParams(location.search);
  var method=q.get('method')||'pfba';
  var cond=q.get('cond')||null;
  var info=document.getElementById('inf-escher-info');
  fetch('infection_flux_by_condition.json?v='+Date.now())
    .then(function(r){if(!r.ok)throw new Error('HTTP '+r.status);return r.json();})
    .then(function(data){
      if(!cond) cond=data.condition_order[0];
      var fx=((data.flux[method]||{})[cond])||{};
      var maxabs=0,nActive=0,k;
      for(k in fx){var a=Math.abs(fx[k]);if(a>1e-6){nActive++;if(a>maxabs)maxabs=a;}}
      var denom=Math.log1p(maxabs)||1;
      function mag(f){var a=Math.abs(f);return a<=1e-6?0:Math.log1p(a)/denom;}
      function paint(){
        var n=0;
        document.querySelectorAll('g.reaction').forEach(function(g){
          var lab=g.querySelector('.reaction-label');if(!lab)return;
          var rxn=lab.textContent.trim().split(/\s+/)[0];
          var f=fx[rxn],active=(f!==undefined&&Math.abs(f)>1e-6),color,size;
          if(active){var t=mag(f);color=interp('#dbe6f2','#0969da',t);size=(2+11*t)+'px';}
          else{color='#e4e8ec';size='1.5px';}
          g.querySelectorAll('path,line,polyline').forEach(function(p){p.style.stroke=color;p.style.strokeWidth=size;});
          lab.style.fill=active?color:'#c4cace';
          n++;
        });
        return n;
      }
      var tries=0,iv=setInterval(function(){
        if(document.querySelectorAll('g.reaction').length>0){clearInterval(iv);paint();}
        else if(++tries>60){clearInterval(iv);}
      },200);
      if(info){
        var lab=(data.labels&&data.labels[cond])||cond;
        var st=(data.stages&&data.stages[cond])||'';
        info.innerHTML='<b>'+esc(method.toUpperCase())+'</b> &middot; '+esc(lab)+
          ' <span style="color:#888">('+esc(st)+')</span><br>'+nActive+
          ' reactions carrying flux &middot; max |v| = '+maxabs.toFixed(2)+' mmol&middot;gDW&#8315;&#185;&middot;h&#8315;&#185;';
      }
    })
    .catch(function(e){if(info)info.textContent='flux data load failed: '+e.message;});
})();
</script>
"""

NEW_LEGEND = """
<div id="atp-safe-legend" style="position:fixed;bottom:20px;right:20px;background:rgba(255,255,255,0.97);border:1px solid #ddd;border-radius:10px;padding:14px 18px;font-family:sans-serif;font-size:12px;z-index:9999;box-shadow:0 3px 12px rgba(0,0,0,0.18);max-width:300px">
  <b style="font-size:13px">FSP237 flux map</b>
  <div id="inf-escher-info" style="margin:8px 0;color:#333;line-height:1.5">loading&hellip;</div>
  <div style="display:flex;align-items:center;gap:8px;margin-top:6px">
    <span style="font-size:11px;color:#888">low</span>
    <span style="flex:1;height:10px;border-radius:3px;background:linear-gradient(90deg,#dbe6f2,#0969da)"></span>
    <span style="font-size:11px;color:#888">high</span>
  </div>
  <i style="color:#999;font-size:11px">edge intensity &amp; width = |flux| (log-scaled); grey = no flux</i>
</div>
"""


def build_flux_json():
    # labels / stages from the unified matrix (pfba aerobic slice)
    long = pd.read_csv(os.path.join(FLUX_RESULTS, "unified_flux_matrix.tsv"),
                       sep="\t", usecols=["condition_id", "O2", "label",
                                          "stage", "method"])
    meta = (long[(long.O2 == "aerobic") & (long.method == "pfba")]
            [["condition_id", "label", "stage"]].drop_duplicates()
            .sort_values("condition_id"))
    cond_order = meta["condition_id"].tolist()
    labels = dict(zip(meta["condition_id"], meta["label"]))
    stages = dict(zip(meta["condition_id"], meta["stage"]))

    flux = {}
    for m in METHODS:
        wide = pd.read_csv(os.path.join(FLUX_RESULTS,
                           f"wide_flux_{m}_aerobic.tsv"), sep="\t")
        wide = wide.set_index("rxn_id")
        per_cond = {}
        for c in cond_order:
            if c not in wide.columns:
                continue
            col = wide[c]
            nz = col[col.abs() > FLUX_EPS]
            per_cond[c] = {rxn: round(float(v), 5) for rxn, v in nz.items()}
        flux[m] = per_cond

    payload = {"methods": METHODS, "condition_order": cond_order,
               "labels": labels, "stages": stages, "flux": flux}
    with open(OUT_JSON, "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    return payload


def build_escher_html():
    html = open(SRC_MAP).read()
    # drop the static single-condition recolor script
    html2 = re.sub(r'<script id="atp-safe-recolor">.*?</script>', "",
                   html, flags=re.S)
    # replace the static legend (through </body>) with dynamic legend + recolor
    repl = NEW_LEGEND + DYNAMIC_RECOLOR + "\n</body>"
    html3, n = re.subn(r'<div id="atp-safe-legend".*?</body>',
                       lambda _m: repl, html2, flags=re.S)
    if n != 1:
        raise RuntimeError(f"legend/body replace matched {n} times (expected 1)")
    # retitle
    html3 = html3.replace("<title>Escher Builder</title>",
                          "<title>FSP237 infection-stage flux map</title>")
    with open(OUT_MAP, "w") as fh:
        fh.write(html3)
    return len(html3)


def main():
    p = build_flux_json()
    size = build_escher_html()
    nnz = sum(len(v) for m in p["flux"].values() for v in m.values())
    print(f"wrote {OUT_JSON}  ({os.path.getsize(OUT_JSON)/1024:.0f} KB, "
          f"{nnz} nonzero flux entries)")
    print(f"wrote {OUT_MAP}  ({size/1024:.0f} KB)")
    # provenance copies
    web = os.path.join(ANALYSIS, "web")
    os.makedirs(web, exist_ok=True)
    for src in (OUT_JSON,):
        with open(src) as fi, open(os.path.join(web, os.path.basename(src)), "w") as fo:
            fo.write(fi.read())


if __name__ == "__main__":
    main()
