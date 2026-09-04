#!/opt/env/modelseed/bin/python
"""Extract experimental kcat (turnover number) per EC from the BRENDA 2025_1
JSON release, organism-tiered (Colletotrichum > other fungi > any). Replaces the
flat EC-class prior used in Stage 2 with real enzyme kinetics.

BRENDA turnover_number values are s^-1. We take the MAXIMUM kcat within the most
specific organism tier that has data (GECKO/AutoPACMEN "catalytic capacity"
convention), plus the count and tier for provenance. EC-hierarchy parents
(3-field, 2-field) are pre-aggregated so a 4-field EC with no direct data can
fall back to its family.

Writes outputs/ec/kcat_brenda_by_ec.tsv:  ec  kcat  tier  n  n_fungal  n_any
"""
import json, re
from pathlib import Path
import numpy as np

ROOT = Path("/home/janakae/fsp237")
BRENDA = Path("/tmp/brenda/brenda_2025_1.json")
OUT = ROOT / "proteomics_integration/outputs/ec/kcat_brenda_by_ec.tsv"

# fungal genera (Ascomycota + Basidiomycota, model-relevant) for the fungal tier
FUNGAL_GENERA = {
    "colletotrichum", "saccharomyces", "candida", "aspergillus", "neurospora",
    "fusarium", "magnaporthe", "pyricularia", "trichoderma", "penicillium",
    "yarrowia", "kluyveromyces", "pichia", "komagataella", "schizosaccharomyces",
    "ustilago", "botrytis", "zymoseptoria", "sclerotinia", "verticillium",
    "cryptococcus", "coprinopsis", "schizophyllum", "pleurotus", "agaricus",
    "rhizopus", "mucor", "mortierella", "phanerochaete", "trametes", "laccaria",
    "malassezia", "debaryomyces", "hansenula", "ogataea", "scheffersomyces",
    "cordyceps", "beauveria", "metarhizium", "chaetomium", "podospora",
    "myceliophthora", "thermothelomyces", "thermomyces", "talaromyces",
    "gibberella", "nectria", "epichloe", "claviceps", "blumeria", "puccinia",
    "melampsora", "sporisorium", "tuber", "yarrowia", "geotrichum",
    "trichophyton", "microsporum", "histoplasma", "blastomyces", "coccidioides",
    "paracoccidioides", "cladosporium", "alternaria", "bipolaris",
    "cochliobolus", "setosphaeria", "leptosphaeria", "phaeosphaeria",
    "dothistroma", "sordaria", "acremonium", "monascus", "wickerhamomyces",
    "lachancea", "torulaspora", "zygosaccharomyces", "eremothecium",
    "ashbya", "cyberlindnera", "starmerella", "brettanomyces", "dekkera",
    "hortaea", "exophiala", "aureobasidium", "sclerotium", "rhizoctonia",
    "thanatephorus", "serpula", "postia", "wolfiporia", "ganoderma",
}
VAL_RE = re.compile(r"^\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)")


def genus(org):
    return org.split()[0].lower() if org else ""


def main():
    d = json.load(open(BRENDA))
    data = d["data"]
    print(f"BRENDA {d.get('version')} — {len(data)} EC entries")

    rows = {}  # ec -> dict
    parent3, parent2 = {}, {}  # hierarchy aggregation (any-organism kcats)
    for ec, e in data.items():
        if not isinstance(e, dict):
            continue
        tns = e.get("turnover_number")
        if not tns:
            continue
        prot = e.get("protein", {}) or {}
        # protein id -> genus
        pgen = {pid: genus(pv.get("organism", "")) for pid, pv in prot.items()}
        fungal, colleto, anyv = [], [], []
        for t in tns:
            m = VAL_RE.match(str(t.get("value", "")))
            if not m:
                continue
            v = float(m.group(1))
            if v <= 0 or v > 1e7:          # skip -999 "not determined" & absurd
                continue
            anyv.append(v)
            gens = {pgen.get(p, "") for p in t.get("proteins", [])}
            if any(g == "colletotrichum" for g in gens):
                colleto.append(v)
            if any(g in FUNGAL_GENERA for g in gens):
                fungal.append(v)
        if not anyv:
            continue
        if colleto:
            kcat, tier = max(colleto), "brenda_colletotrichum"
        elif fungal:
            kcat, tier = max(fungal), "brenda_fungal"
        else:
            kcat, tier = np.percentile(anyv, 90), "brenda_any"   # 90th pct: capacity, outlier-robust
        rows[ec] = dict(ec=ec, kcat=round(float(kcat), 4), tier=tier,
                        n=len(anyv), n_fungal=len(fungal), n_any=len(anyv))
        # feed hierarchy parents with all any-organism values
        parts = ec.split(".")
        if len(parts) >= 3:
            parent3.setdefault(".".join(parts[:3]), []).extend(anyv)
        if len(parts) >= 2:
            parent2.setdefault(".".join(parts[:2]), []).extend(anyv)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        fh.write("ec\tkcat\ttier\tn\tn_fungal\tn_any\n")
        for ec in sorted(rows):
            r = rows[ec]
            fh.write(f"{r['ec']}\t{r['kcat']}\t{r['tier']}\t{r['n']}\t{r['n_fungal']}\t{r['n_any']}\n")

    # hierarchy parent medians (for EC with no direct BRENDA data)
    p3 = ROOT / "proteomics_integration/outputs/ec/kcat_brenda_parent3.tsv"
    with open(p3, "w") as fh:
        fh.write("ec3\tkcat_median\tn\n")
        for k in sorted(parent3):
            vs = parent3[k]
            fh.write(f"{k}\t{round(float(np.median(vs)),4)}\t{len(vs)}\n")

    tiers = {}
    for r in rows.values():
        tiers[r["tier"]] = tiers.get(r["tier"], 0) + 1
    print(f"EC with kcat: {len(rows)}  | tiers: {tiers}")
    print(f"parent-3 families: {len(parent3)}")
    print(f"wrote {OUT.name}, {p3.name}")


if __name__ == "__main__":
    main()
