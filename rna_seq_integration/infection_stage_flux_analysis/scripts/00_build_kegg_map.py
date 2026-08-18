"""00 - Build the ModelSEED -> KEGG map for the hybrid pathway taxonomy.

Primary pathway grouping remains the validated 23-bucket classifier
(outputs/pathway_assignment.tsv). This step ATTACHES, where derivable, real KEGG
reaction IDs, KEGG pathway IDs/names, and EC numbers, keyed on the model rxn_id.

Source (read-only): the ModelSEEDDatabase clone. We snapshot the derived map into
data/rxn_kegg_map.tsv so this analysis is self-contained if that path moves.

Coverage is expected to be partial (~44% KEGG reaction, ~41% KEGG pathway): the
1622 model reactions include 337 custom ids (EX_/bio/frxn/BiGG-style) with no
ModelSEED base and thus no KEGG mapping -- those keep the 23-bucket label only.
"""
import os
import re
import pandas as pd
from _common import OUTPUTS, DATA, MODELSEED_DB, seed_base, ensure_dirs

ALIASES = os.path.join(MODELSEED_DB, "Aliases", "Unique_ModelSEED_Reaction_Aliases.txt")
REACTIONS = os.path.join(MODELSEED_DB, "reactions.tsv")

KEGG_PW_RE = re.compile(r"rn(\d+)\s*\(([^)]+)\)")


def parse_kegg_pathways(pathways_field: str):
    """From the reactions.tsv 'pathways' field, pull only the KEGG segment."""
    if not isinstance(pathways_field, str) or not pathways_field or pathways_field == "null":
        return [], []
    # field looks like:  MetaCyc: ...|KEGG: rn00230 (Purine metabolism); rn...
    kegg_seg = ""
    for seg in pathways_field.split("|"):
        seg = seg.strip()
        if seg.startswith("KEGG:"):
            kegg_seg = seg[len("KEGG:"):]
            break
    ids, names = [], []
    for m in KEGG_PW_RE.finditer(kegg_seg):
        ids.append("rn" + m.group(1))
        names.append(m.group(2).strip())
    return ids, names


def main():
    ensure_dirs()

    # all model reactions (authoritative id list + primary pathway bucket)
    pa = pd.read_csv(os.path.join(OUTPUTS, "pathway_assignment.tsv"), sep="\t")
    pa["seed_base"] = pa["rxn_id"].map(seed_base)

    # 1) KEGG reaction ids from the alias table (Source == KEGG)
    al = pd.read_csv(ALIASES, sep="\t", names=["seed", "ext", "source"], header=0)
    kegg_rxn = (al[al["source"] == "KEGG"]
                .groupby("seed")["ext"]
                .apply(lambda s: ";".join(sorted(set(s))))
                .to_dict())

    # 2) KEGG pathways + EC from reactions.tsv
    rx = pd.read_csv(REACTIONS, sep="\t", low_memory=False)
    kegg_pw_id, kegg_pw_name, ec_map = {}, {}, {}
    for _, r in rx.iterrows():
        rid = r["id"]
        ids, names = parse_kegg_pathways(r.get("pathways", ""))
        if ids:
            kegg_pw_id[rid] = ";".join(ids)
            kegg_pw_name[rid] = ";".join(dict.fromkeys(names))  # dedupe, keep order
        ec = r.get("ec_numbers", "")
        if isinstance(ec, str) and ec and ec != "null":
            ec_map[rid] = ec

    def g(d, base):
        return d.get(base, "") if isinstance(base, str) else ""

    out = pd.DataFrame({
        "rxn_id": pa["rxn_id"],
        "seed_base": pa["seed_base"].fillna(""),
        "pathway_bucket": pa["pathway"],
        "kegg_reaction": pa["seed_base"].map(lambda b: g(kegg_rxn, b)),
        "kegg_pathway_id": pa["seed_base"].map(lambda b: g(kegg_pw_id, b)),
        "kegg_pathway_name": pa["seed_base"].map(lambda b: g(kegg_pw_name, b)),
        "ec_number": pa["seed_base"].map(lambda b: g(ec_map, b)),
    })

    out_path = os.path.join(DATA, "rxn_kegg_map.tsv")
    out.to_csv(out_path, sep="\t", index=False)

    n = len(out)
    n_base = int((out["seed_base"] != "").sum())
    n_kr = int((out["kegg_reaction"] != "").sum())
    n_kp = int((out["kegg_pathway_id"] != "").sum())
    n_ec = int((out["ec_number"] != "").sum())
    n_distinct_pw = len({p for cell in out["kegg_pathway_id"] if cell
                         for p in cell.split(";")})
    print(f"wrote {out_path}")
    print(f"  reactions ................ {n}")
    print(f"  with ModelSEED base ...... {n_base} ({100*n_base/n:.1f}%)")
    print(f"  with KEGG reaction id .... {n_kr} ({100*n_kr/n:.1f}%)")
    print(f"  with KEGG pathway ........ {n_kp} ({100*n_kp/n:.1f}%)  "
          f"[{n_distinct_pw} distinct KEGG pathways]")
    print(f"  with EC number ........... {n_ec} ({100*n_ec/n:.1f}%)")


if __name__ == "__main__":
    main()
