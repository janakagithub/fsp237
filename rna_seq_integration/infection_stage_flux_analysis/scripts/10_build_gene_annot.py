"""10 - Gene annotation, pathogenicity flags, and candidate-effector layer.

Answers the user question "are any of the genes behind active-flux reactions
effectors related to pathogenesis?" honestly, with three layers:

  1. by_reaction        : reaction -> its GPR (metabolic) genes, each with a real
                          function name (from the proteome), S1 expression + bin.
  2. virulence_metabolic: the subset of GPR (metabolic) genes whose enzyme family
                          is a known virulence/pathogenicity factor (pectinases,
                          cutinase, glyoxylate cycle, trehalose, cell-wall synthases,
                          ROS detox, melanin, P450). These are metabolic virulence
                          factors, NOT canonical effectors.
  3. candidate_effectors: highly-expressed genes that are NOT in the metabolic model
                          (non-metabolic) and look secreted/effector-like. This is the
                          honest home for "effectors": canonical fungal effectors are
                          small secreted proteins that (by construction) sit outside a
                          metabolic GEM's gene set, so they surface here, not as GPR
                          genes.

Data sources (self-contained under fsp237/):
  - proteome  : BLAST protein DB gpr-update/blast_db/fsp237  (titles carry
                `gene_NNNN len=.. func=<description>`), extracted via blastdbcmd.
  - expression: expression-data/S1_normalized_expression.xlsx (13,047 genes).
  - GPR/model : outputs/reaction_expression.tsv (gpr column -> model gene set).

Secretome prediction: SignalP 6.0 / EffectorP 3.0 are NOT runnable in this
environment (no EMBOSS pepstats, no torch, no academic license). We therefore use
a transparent, reproducible sequence-feature HEURISTIC (N-terminal hydrophobic
signal-peptide core + small + cysteine-rich + non-metabolic function), clearly
labelled as such in the payload/README. All other layers are unaffected.
"""
import os
import re
import sys
import json
import subprocess
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from _common import ANALYSIS  # noqa: E402

FSP = "/home/janakae/fsp237"
PROT_DB = os.path.join(FSP, "gpr-update", "blast_db", "fsp237")
EXPR_XLSX = os.path.join(FSP, "expression-data", "S1_normalized_expression.xlsx")
REXP = os.path.join(FSP, "rna_seq_integration", "outputs", "reaction_expression.tsv")
DATA = os.path.join(ANALYSIS, "data")
WEB = os.path.join(ANALYSIS, "web")
FAA = os.path.join(DATA, "fsp237_proteome.faa")
ANNOT_TSV = os.path.join(DATA, "gene_annotation.tsv")
SEC_TSV = os.path.join(DATA, "secretome_predictions.tsv")
OUT_JSON = "/home/janakae/fsp237/atp-safe/infection_gene_annot.json"

GENE_RE = re.compile(r"(?:gene_\d+|CH63R_\S+)")
FUNC_RE = re.compile(r"func=(.*)$")

# Kyte-Doolittle hydropathy (signal-peptide h-region detection)
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5,
      'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
      'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
      'Y': -1.3, 'V': 4.2}

# Metabolic virulence / pathogenicity enzyme families (keyword match on func=).
# Sourced from the M3 deep-research literature on Colletotrichum pathogenesis.
VIRU_FAMILIES = [
    ("Pectin/cell-wall degradation",
     ["pectate lyase", "pectin lyase", "polygalacturonase", "pectinesterase",
      "pectin methylesterase", "rhamnogalacturon", "galacturonase",
      "arabinofuranosidase", "xylanase", "cellulase", "glucanase",
      "cellobiohydrolase", "endoglucanase"]),
    ("Cutinase / lipase (surface breach)", ["cutinase", "lipase", "esterase"]),
    ("Glyoxylate cycle", ["isocitrate lyase", "malate synthase"]),
    ("Trehalose metabolism (appressorium turgor)",
     ["trehalose-6-phosphate synthase", "trehalose synthase", "tps1",
      "trehalose-phosphatase", "trehalose phosphatase", "trehalase"]),
    ("Fungal cell-wall synthesis (chitin/glucan)",
     ["chitin synthase", "glucan synthase", "1,3-beta-glucan",
      "beta-1,3-glucan", "chitinase"]),
    ("Peroxisomal beta-oxidation",
     ["acyl-coa oxidase", "enoyl-coa hydratase", "3-ketoacyl",
      "3-hydroxyacyl-coa", "carnitine"]),
    ("ROS detoxification", ["catalase", "superoxide dismutase", "peroxidase",
                            "peroxiredoxin", "glutathione peroxidase"]),
    ("Melanin / pigment (appressorium)",
     ["laccase", "tyrosinase", "polyketide synthase", "scytalone",
      "1,3,8-trihydroxynaphthalene", "melanin"]),
    ("Cytochrome P450 (detox/secondary metab)", ["cytochrome p450"]),
]

# Effector/secretion function keywords (non-enzymatic small-secreted signatures).
EFFECTOR_KW = ["effector", "secreted", "hydrophobin", "lysm", "cerato-platanin",
               "necrosis", "necrosis-inducing", "snodprot", "cfem", "avr",
               "elicitin", "npp1", "ecp", "small secreted"]

# Function terms that mark a protein as clearly metabolic/housekeeping (exclude
# from effector candidacy even if small).
METABOLIC_KW = ["synthase", "synthetase", "dehydrogenase", "reductase",
                "oxidase", "transferase", "kinase", "phosphatase", "hydrolase",
                "isomerase", "lyase", "ligase", "carboxylase", "decarboxylase",
                "aminotransferase", "hydratase", "mutase", "polymerase",
                "transporter", "permease", "atpase", "ribosomal", "elongation factor",
                "trna", "rrna", "helicase", "topoisomerase", "aldolase",
                "epimerase", "racemase", "thioesterase", "acyltransferase"]


def extract_proteome():
    """blastdbcmd -> FASTA with real gene_NNNN headers (len=/func= in title)."""
    os.makedirs(DATA, exist_ok=True)
    with open(FAA, "w") as fh:
        subprocess.run(["blastdbcmd", "-db", PROT_DB, "-entry", "all",
                        "-outfmt", "%f"], check=True, stdout=fh)
    # parse fasta (id, description, sequence) without a hard Biopython dep
    seqs, func, desc = {}, {}, {}
    gid, buf = None, []
    with open(FAA) as fh:
        for line in fh:
            if line.startswith(">"):
                if gid is not None:
                    seqs[gid] = "".join(buf)
                header = line[1:].strip()
                gid = header.split()[0]
                m = FUNC_RE.search(header)
                func[gid] = m.group(1).strip() if m else ""
                desc[gid] = header
                buf = []
            else:
                buf.append(line.strip())
        if gid is not None:
            seqs[gid] = "".join(buf)
    return seqs, func


def signal_peptide(seq):
    """Conservative N-terminal signal-peptide heuristic (SignalP-lite):
    a hydrophobic h-region core in the first ~25 residues."""
    s = seq[:30]
    if len(s) < 12:
        return False
    best = -9.9
    for start in range(1, min(16, len(s) - 8)):
        win = s[start:start + 9]
        vals = [KD.get(a, 0.0) for a in win]
        best = max(best, sum(vals) / len(vals))
    # require a strong hydrophobic core and a not-fully-charged n-region
    nreg = s[:5]
    n_charged = sum(1 for a in nreg if a in "DEKR")
    return best >= 1.6 and n_charged <= 2


def secretome(seqs, func):
    rows = {}
    for g, seq in seqs.items():
        L = len(seq)
        cys = seq.count("C")
        cys_frac = cys / L if L else 0.0
        sig = signal_peptide(seq)
        f = (func.get(g, "") or "").lower()
        kw_eff = any(k in f for k in EFFECTOR_KW)
        is_metab = any(k in f for k in METABOLIC_KW)
        eff = ((sig and L <= 300 and (cys_frac >= 0.02 or cys >= 4) and not is_metab)
               or kw_eff)
        rows[g] = {
            "length": L, "cys": cys, "cys_frac": round(cys_frac, 4),
            "sig_peptide": bool(sig), "effector_like": bool(eff),
            "kw_effector": bool(kw_eff),
            "localization": "secreted (pred.)" if sig or kw_eff else "other",
        }
    return rows


def viru_family(func_str):
    f = (func_str or "").lower()
    for fam, kws in VIRU_FAMILIES:
        if any(k in f for k in kws):
            return fam
    return ""


def main():
    os.makedirs(DATA, exist_ok=True)
    os.makedirs(WEB, exist_ok=True)

    print("extracting proteome via blastdbcmd ...", flush=True)
    seqs, func = extract_proteome()
    print(f"  proteome: {len(seqs)} proteins", flush=True)

    # expression
    xl = pd.read_excel(EXPR_XLSX, sheet_name="Normalized log2 TPM").set_index("Gene ID")
    mean_by_gene = xl["Mean log2(TPM+1)"].astype(float).to_dict()
    mean_by_gene = {str(k): v for k, v in mean_by_gene.items()}
    all_expr_genes = list(mean_by_gene.keys())
    # full-transcriptome percentile rank
    ser = pd.Series(mean_by_gene)
    pct = ser.rank(pct=True).to_dict()

    # GPR / model gene set + reaction -> genes
    rexp = pd.read_csv(REXP, sep="\t", comment="#")
    rxn_genes = {}
    model_genes = set()
    for _, r in rexp.iterrows():
        gpr = str(r.get("gpr", "") or "")
        genes = sorted(set(GENE_RE.findall(gpr)))
        if genes:
            rxn_genes[r["rxn_id"]] = genes
            model_genes.update(genes)
    print(f"  model genes (in a GPR): {len(model_genes)}", flush=True)

    # bin cutoffs on model-gene distribution (matches src/gpr_expression.py)
    ms = pd.Series({g: mean_by_gene.get(g, np.nan) for g in model_genes}).dropna()
    q25, q75 = float(ms.quantile(0.25)), float(ms.quantile(0.75))

    def gbin(g):
        v = mean_by_gene.get(g)
        if v is None or (isinstance(v, float) and np.isnan(v)) or v <= 0:
            return "absent"
        if v >= q75:
            return "hi"
        if v >= q25:
            return "med"
        return "lo"

    sec = secretome(seqs, func)

    def gene_row(g, in_model):
        v = mean_by_gene.get(g)
        s = sec.get(g, {})
        return {
            "gene": g,
            "func": func.get(g, "") or "(no annotation)",
            "log2tpm": None if v is None or (isinstance(v, float) and np.isnan(v)) else round(float(v), 3),
            "bin": gbin(g),
            "in_model": in_model,
            "secreted": bool(s.get("sig_peptide") or s.get("kw_effector")),
            "effector": bool(s.get("effector_like")),
            "viru": viru_family(func.get(g, "")),
            "len": s.get("length"),
        }

    # layer 1: by_reaction
    by_reaction = {rxn: [gene_row(g, True) for g in genes]
                   for rxn, genes in rxn_genes.items()}

    # layer 2: metabolic virulence factors (GPR genes in a virulence family)
    gene_to_rxns = {}
    for rxn, genes in rxn_genes.items():
        for g in genes:
            gene_to_rxns.setdefault(g, []).append(rxn)
    virulence_metabolic = []
    for g in sorted(model_genes):
        fam = viru_family(func.get(g, ""))
        if not fam:
            continue
        v = mean_by_gene.get(g)
        virulence_metabolic.append({
            "gene": g, "func": func.get(g, "") or "(no annotation)",
            "family": fam,
            "log2tpm": None if v is None or (isinstance(v, float) and np.isnan(v)) else round(float(v), 3),
            "bin": gbin(g),
            "n_rxn": len(gene_to_rxns.get(g, [])),
            "rxns": sorted(gene_to_rxns.get(g, []))[:12],
        })
    virulence_metabolic.sort(key=lambda d: (d["log2tpm"] is None, -(d["log2tpm"] or 0)))

    # layer 3: candidate effectors = NON-metabolic (not in model) genes that are
    # highly expressed and look secreted/effector-like.
    non_model = [g for g in all_expr_genes if g not in model_genes]
    cand = []
    for g in non_model:
        s = sec.get(g, {})
        v = mean_by_gene.get(g)
        if v is None or (isinstance(v, float) and np.isnan(v)):
            continue
        p = pct.get(g, 0.0)
        high = p >= 0.90
        secreted = bool(s.get("sig_peptide") or s.get("kw_effector"))
        effector = bool(s.get("effector_like"))
        fl = (func.get(g, "") or "").lower()
        is_metab = any(k in fl for k in METABOLIC_KW)
        # keep genes that are effector-like, OR (secreted AND highly expressed and
        # not a clearly-metabolic/housekeeping protein), OR an effector keyword hit.
        if effector or (secreted and high and not is_metab) or s.get("kw_effector"):
            cand.append({
                "gene": g, "func": func.get(g, "") or "(no annotation)",
                "log2tpm": round(float(v), 3), "pctile": round(float(p), 3),
                "len": s.get("length"),
                "secreted": secreted, "effector": effector, "high_expr": bool(high),
            })
    cand.sort(key=lambda d: -d["log2tpm"])
    cand = cand[:250]

    # annotation + secretome TSV provenance
    with open(ANNOT_TSV, "w") as fh:
        fh.write("gene_id\tlength\tfunc\n")
        for g in sorted(seqs):
            fh.write(f"{g}\t{len(seqs[g])}\t{func.get(g,'')}\n")
    with open(SEC_TSV, "w") as fh:
        fh.write("gene_id\tlength\tcys\tcys_frac\tsig_peptide\teffector_like\t"
                 "kw_effector\tlocalization\n")
        for g in sorted(sec):
            s = sec[g]
            fh.write(f"{g}\t{s['length']}\t{s['cys']}\t{s['cys_frac']}\t"
                     f"{int(s['sig_peptide'])}\t{int(s['effector_like'])}\t"
                     f"{int(s['kw_effector'])}\t{s['localization']}\n")

    n_sig = sum(1 for s in sec.values() if s["sig_peptide"])
    n_eff = sum(1 for s in sec.values() if s["effector_like"])

    payload = {
        "meta": {
            "note": ("Reaction->gene drill-down and pathogenicity cross-reference. "
                     "The model's GPR genes are metabolic ENZYMES; canonical fungal "
                     "effectors are small secreted proteins that, by construction, sit "
                     "OUTSIDE a metabolic model's gene set. So 'metabolic virulence "
                     "factors' (enzyme families with known roles in pathogenesis) are "
                     "reported separately from 'candidate effectors' (highly-expressed "
                     "NON-metabolic secreted/effector-like genes)."),
            "secretome_method": ("Sequence-feature HEURISTIC (N-terminal hydrophobic "
                                 "signal-peptide core + small + cysteine-rich + "
                                 "non-metabolic function). SignalP 6.0 / EffectorP 3.0 "
                                 "were not runnable in this environment (no EMBOSS "
                                 "pepstats, no torch, no academic license); treat "
                                 "secreted/effector flags as indicative, not definitive."),
            "bin_cutoffs": {"hi_ge": round(q75, 3), "lo_ge": round(q25, 3)},
            "n_proteins": len(seqs),
            "n_expr_genes": len(all_expr_genes),
            "n_model_genes": len(model_genes),
            "n_sig_peptide": n_sig,
            "n_effector_like": n_eff,
            "n_candidate_effectors": len(cand),
            "n_virulence_metabolic": len(virulence_metabolic),
        },
        "by_reaction": by_reaction,
        "virulence_metabolic": virulence_metabolic,
        "candidate_effectors": cand,
    }

    with open(OUT_JSON, "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    with open(os.path.join(WEB, "infection_gene_annot.json"), "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))

    print(f"wrote {OUT_JSON} ({os.path.getsize(OUT_JSON)/1024:.0f} KB)")
    print(f"  reactions with GPR : {len(by_reaction)}")
    print(f"  virulence_metabolic: {len(virulence_metabolic)} "
          f"(top: {virulence_metabolic[0]['func'][:50] if virulence_metabolic else '-'})")
    print(f"  candidate_effectors: {len(cand)} (non-metabolic, high-expr/secreted)")
    print(f"  secretome heuristic: {n_sig} signal-peptide, {n_eff} effector-like "
          f"of {len(seqs)} proteins")


if __name__ == "__main__":
    main()
