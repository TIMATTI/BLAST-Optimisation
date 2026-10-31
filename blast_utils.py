# blast_utils.py
from pathlib import Path
import subprocess
import shutil
from decimal import Decimal
from typing import List, Tuple, Dict
import sys

# безопасный импорт pandas
try:
    import pandas as pd
except Exception:
    pd = None

# OUTFMT: используем slen и bitscore (надежные поля)
OUTFMT_FIELDS = [
    "qseqid","sseqid","pident","length","mismatch","gapopen",
    "qstart","qend","sstart","send","slen","evalue","bitscore","score"
]
OUTFMT = "6 " + ",".join(OUTFMT_FIELDS)

# ----------------- Tools check -----------------
def check_tools():
    for t in ("blastn", "makeblastdb"):
        if shutil.which(t) is None:
            raise RuntimeError(f"Required program '{t}' not found in PATH.")

# ----------------- FASTA helpers -----------------
def fasta_wrap(seq: str, width: int = 80) -> str:
    seq = "".join(seq.split())
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def write_fasta_from_seq(seq: str, out_path: Path, header: str = "seq") -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write(f">{header}\n{fasta_wrap(seq)}\n")
    return out_path

def makeblastdb(fasta_path: Path, db_out_prefix: Path) -> Path:
    cmd = ["makeblastdb", "-in", str(fasta_path), "-dbtype", "nucl", "-out", str(db_out_prefix)]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return db_out_prefix

# ----------------- Run BLAST -----------------
def run_blast(db_prefix_or_subject: str, query_path: str, out_path: Path,
              use_db: bool=True, threads: int=1, evalue: float=1e-3, max_target_seqs: int=100) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["blastn", "-query", str(query_path)]
    if use_db:
        cmd += ["-db", db_prefix_or_subject]
    else:
        cmd += ["-subject", db_prefix_or_subject]
    cmd += [
        "-outfmt", OUTFMT,
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(threads),
        "-out", str(out_path)
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return out_path

# ----------------- Load results -----------------
def load_results_tsv(tsv_path: Path):
    """
    Robust loader: auto-detects number of columns and chooses matching field names.
    Returns pandas.DataFrame (if pandas available) or list of dicts.
    """
    OLD12 = ["qseqid","sseqid","pident","length","mismatch","gapopen",
             "qstart","qend","sstart","send","evalue","score"]
    OLD13 = ["qseqid","sseqid","pident","length","mismatch","gapopen",
             "qstart","qend","sstart","send","evalue","bitscore","score"]
    # OUTFMT_FIELDS from module (14 fields) may be present as well

    lines = []
    with tsv_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line:  # skip empty lines
                lines.append(line)

    if not lines:
        return pd.DataFrame() if pd else []

    # detect columns from first non-empty line
    first_parts = [p.strip() for p in lines[0].split("\t")]
    ncols = len(first_parts)

    if ncols == len(OUTFMT_FIELDS):
        fields = OUTFMT_FIELDS
    elif ncols == 12:
        fields = OLD12
    elif ncols == 13:
        fields = OLD13
    else:
        # fallback: assume first N fields of OUTFMT_FIELDS
        fields = OUTFMT_FIELDS[:ncols]

    rows = []
    for line in lines:
        parts = [p.strip() for p in line.split("\t")]
        # pad or truncate to match fields length
        if len(parts) < len(fields):
            parts += [""] * (len(fields) - len(parts))
        elif len(parts) > len(fields):
            # keep extras as extra_0...
            pass
        row = {k: (parts[i] if i < len(parts) else "") for i, k in enumerate(fields)}
        # store extras if any
        if len(parts) > len(fields):
            for j, extra in enumerate(parts[len(fields):]):
                row[f"extra_{j}"] = extra
        rows.append(row)

    if pd:
        df = pd.DataFrame(rows)
        for col in ("pident","evalue","bitscore","score","slen","length","qstart","qend","sstart","send"):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        # debug: print detected mapping once (optional)
        # print("DEBUG: detected TSV fields:", fields, file=sys.stderr)
        return df
    return rows

