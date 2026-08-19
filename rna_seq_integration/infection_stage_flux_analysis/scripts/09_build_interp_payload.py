"""09 - Build the M3 "Biological Interpretation" payload.

Reads the deep-research markdown files under deep_research/pathways/ (one per
priority pathway, authored by the M3 research agents), extracts each file's
closing ```json summary block, renders the prose (everything before that block)
to HTML with mistune, and emits atp-safe/infection_interp.json for the new
"Biological Interpretation" sub-tab.

No analysis is recomputed here; this only packages the narrative. Run after the
research markdown is in place.
"""
import os
import re
import json
import html
import mistune
from _common import ANALYSIS

PATH_DIR = os.path.join(ANALYSIS, "deep_research", "pathways")
OUT_JSON = "/home/janakae/fsp237/atp-safe/infection_interp.json"

# display order + titles + how each pathway entered the priority list
REGISTRY = [
    ("plant_cell_wall_degradation",
     "Plant cell-wall degradation (pectin & pentoses)", "stage-differential",
     "Most stage-differential — active only in necrotrophic-stage media."),
    ("peroxisomal_beta_oxidation",
     "Peroxisomal β-oxidation", "both",
     "Stage-differential (pre-infection only) and 2nd most method-discordant."),
    ("storage_glycogen_trehalose",
     "Storage carbohydrates (glycogen / trehalose)", "both",
     "Mildly pre-infection-weighted; moderately method-discordant."),
    ("cell_wall_polysaccharide",
     "Fungal cell-wall polysaccharide", "stage-differential",
     "Slightly biotrophic-weighted cell-wall remodeling."),
    ("central_carbon_metabolism",
     "Central carbon (glycolysis / TCA / PPP)", "method-discordant",
     "Uniformly active across stages; where integration methods most disagree."),
    ("gam_maintenance",
     "GAM / maintenance energy", "method-discordant",
     "The single most method-discordant bucket — a non-gene-associated artifact."),
]

JSON_BLOCK = re.compile(r"```json\s*(\{.*?\})\s*```", re.S)


def parse_file(slug):
    fp = os.path.join(PATH_DIR, slug + ".md")
    if not os.path.exists(fp):
        return None
    raw = open(fp).read()
    m = JSON_BLOCK.search(raw)
    summary = {}
    if m:
        try:
            summary = json.loads(m.group(1))
        except json.JSONDecodeError:
            summary = {}
    prose = JSON_BLOCK.sub("", raw).strip()
    html = mistune.html(prose)
    return {"summary": summary, "html": html}


def main():
    pathways = []
    missing = []
    for slug, title, category, why in REGISTRY:
        parsed = parse_file(slug)
        if parsed is None:
            missing.append(slug)
            continue
        s = parsed["summary"]
        # agents sometimes HTML-escape (&amp;, &rho;) inside JSON strings; decode
        # to raw text so the client-side esc() escapes exactly once.
        refs = []
        for r in s.get("top_refs", []):
            refs.append({"cite": html.unescape(str(r.get("cite", ""))),
                         "url": r.get("url", "")})
        pathways.append({
            "slug": slug,
            "title": title,
            "category": category,
            "why": why,
            "one_liner": html.unescape(str(s.get("one_liner", ""))),
            "confidence": html.unescape(str(s.get("confidence", ""))),
            "top_refs": refs,
            "html": parsed["html"],
        })

    payload = {
        "meta": {
            "note": ("Deep-research biological interpretation (M3) of the most "
                     "stage-differential and most method-discordant pathways from "
                     "M1/M2. Evidence is tiered (Tier 1 Colletotrichum-specific, "
                     "Tier 2 related phytopathogens, Tier 3 general fungal/"
                     "biochemical). Stage-differential signals are flux/medium-"
                     "driven: a single static S1 transcriptome is applied to every "
                     "medium, so expression does not vary by stage."),
            "n_pathways": len(pathways),
        },
        "pathways": pathways,
    }

    with open(OUT_JSON, "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    web = os.path.join(ANALYSIS, "web")
    os.makedirs(web, exist_ok=True)
    with open(os.path.join(web, "infection_interp.json"), "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))

    print(f"wrote {OUT_JSON} ({os.path.getsize(OUT_JSON)/1024:.0f} KB, "
          f"{len(pathways)} pathways)")
    if missing:
        print("  MISSING markdown for:", ", ".join(missing))
    for p in pathways:
        print(f"  - {p['slug']:32s} conf={p['confidence']:8s} "
              f"refs={len(p['top_refs'])}")


if __name__ == "__main__":
    main()
