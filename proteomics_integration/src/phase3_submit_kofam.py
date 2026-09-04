#!/usr/bin/env python3
"""Phase 3 (kickoff) — submit full FSP237 proteome to kofam_scan on poplar via Celery.

Annotates every protein with KEGG KO (kofamscan 'mapper' output). KO -> EC then feeds
the kcat parameterization (BRENDA/SABIO EC hierarchy + DLKcat gap-fill) for the
enzyme-constrained (sMOMENT/GECKO) models. One full-proteome run is reused across
kcat, effector, and orphan-resolution analyses.

Broker: redis://bioseed_redis:6379/10 (reachable from seed.jupyter). Progress:
http://poplar.cels.anl.gov:5555/  (Flower).
"""
import os, json, datetime
from pathlib import Path
from Bio import SeqIO
from celery import Celery

ROOT = Path("/home/janakae/fsp237")
FASTA = ROOT / "rna_seq_integration/infection_stage_flux_analysis/data/fsp237_proteome.faa"
OUT = ROOT / "proteomics_integration" / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

JOB_NAME = "fsp237_fullproteome_kofam_20260903"   # unique -> /scratch/bioseed_tools/kofam_scan/run/<JOB_NAME>
KO_DB_DATE = "2025-11-03"
THREADS = 40

# ---- build {protein_id: sequence} ----
proteins = {r.id: str(r.seq) for r in SeqIO.parse(str(FASTA), "fasta")}
print(f"[kofam] proteins to submit: {len(proteins)}")

# ---- celery client (default broker matches poplar's bioseed setup) ----
BROKER = os.getenv("CELERY_BROKER_URL", "redis://bioseed_redis:6379/10")
BACKEND = os.getenv("CELERY_RESULT_BACKEND", "redis://bioseed_redis:6379/10")
bioseed_exec = Celery("client", broker=BROKER, backend=BACKEND)

task = bioseed_exec.send_task(
    "kofam_scan.exec_annotation",
    args=[THREADS, JOB_NAME, KO_DB_DATE, proteins, "mapper"],
    queue="kofam_scan",
)
print(f"[kofam] submitted task id: {task.id}")

meta = dict(task_id=task.id, job_name=JOB_NAME, ko_db_date=KO_DB_DATE,
            threads=THREADS, n_proteins=len(proteins), broker=BROKER,
            output_hint=f"/scratch/bioseed_tools/kofam_scan/run/{JOB_NAME}",
            flower="http://poplar.cels.anl.gov:5555/",
            submitted=datetime.datetime.now().isoformat(timespec="seconds"))
(OUT / "kofam_job.json").write_text(json.dumps(meta, indent=2))
print("[kofam] wrote", OUT / "kofam_job.json")
print(json.dumps(meta, indent=2))
