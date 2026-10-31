#!/usr/bin/env python3
"""
BLAST.py — wrapper для локального BLAST+ (blastn) с подсчётами как у NCBI BLAST.
Подходит для Windows, выводит Accession Length из FASTA и корректные e-values.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Tuple, Dict
from decimal import Decimal, getcontext

getcontext().prec = 12  # точность для суммирования score

try:
    import pandas as pd
except ImportError:
    pd = None

# ----------------- Outfmt -----------------
OUTFMT_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "score"
]
OUTFMT = "6 " + ",".join(OUTFMT_FIELDS)

# ----------------- Tools check -----------------
def check_tools():
    """Проверяет наличие blastn и makeblastdb в PATH."""
    for t in ("blastn", "makeblastdb"):
        if shutil.which(t) is None:
            print(f"Error: required program '{t}' not found in PATH. Install BLAST+.", file=sys.stderr)
            sys.exit(2)

# ----------------- FASTA helpers -----------------
def fasta_wrap(seq: str, width: int = 80) -> str:
    """Разбивает последовательность на строки фиксированной ширины."""
    seq = "".join(seq.split())
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def write_fasta_from_seq(seq: str, out_path: Path, header: str = "seq") -> Path:
    """Записывает последовательность в FASTA файл."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write(f">{header}\n{fasta_wrap(seq)}\n")
    return out_path

def makeblastdb(fasta_path: Path, db_out_prefix: Path) -> Path:
    """Создаёт локальную BLAST-базу из FASTA."""
    cmd = ["makeblastdb", "-in", str(fasta_path), "-dbtype", "nucl", "-out", str(db_out_prefix)]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return db_out_prefix

# ----------------- Prepare query and subjects -----------------
def prepare_query_file(query_arg: str, workdir: Path) -> Tuple[str, int]:
    """
    Подготавливает query: возвращает путь к файлу FASTA и длину последовательности.
    Если query_arg — это файл, берём его, иначе создаём временный FASTA.
    """
    p = Path(query_arg)
    if p.exists():
        text = p.read_text()
        seq = "".join(text.split("\n")[1:]) if text.startswith(">") else text
        return str(p), len(seq)
    qf = workdir / "query.fasta"
    write_fasta_from_seq(query_arg, qf, header="query")
    return str(qf), len("".join(query_arg.split()))

def prepare_combined_subjects(subject_args: List[str], workdir: Path) -> Tuple[str, List[Path], Dict[str,int]]:
    """
    Подготавливает объединённый FASTA файл для субъектов, создаёт BLAST-базу.
    Возвращает префикс базы, список созданных файлов и словарь Accession->Length.
    """
    created = []
    seqs = []
    acc_len_dict = {}
    for i, subj in enumerate(subject_args):
        p = Path(subj)
        if p.exists():
            text = p.read_text()
            if text.lstrip().startswith(">"):
                header = text.split("\n", 1)[0][1:]
                seq = "".join(text.split("\n")[1:])
            else:
                header = f"subject_file_{i}"
                seq = text
        else:
            header = f"subject_{i}"
            seq = subj
        seqs.append(f">{header}\n{fasta_wrap(seq)}")
        acc_len_dict[header] = len(seq)

    combined_fasta = workdir / "combined_subjects.fasta"
    with combined_fasta.open("w") as fh:
        for s in seqs:
            fh.write(s.strip() + "\n")
    created.append(combined_fasta)

    db_prefix = workdir / "subjects_db"
    makeblastdb(combined_fasta, db_prefix)
    created.append(db_prefix)
    return str(db_prefix), created, acc_len_dict

# ----------------- Run BLAST -----------------
def run_blast(db_prefix_or_subject: str, query_path: str, out_path: Path,
              use_db: bool=True, threads: int=1, evalue: float=1e-3, max_target_seqs: int=100) -> Path:
    """
    Запускает blastn на локальной базе или по субъекту.
    """
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
    """Загружает TSV результат BLAST в DataFrame или список словарей."""
    rows = []
    with tsv_path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(OUTFMT_FIELDS):
                parts += [""] * (len(OUTFMT_FIELDS) - len(parts))
            rows.append({k: v for k, v in zip(OUTFMT_FIELDS, parts)})
    if pd:
        df = pd.DataFrame(rows)
        for col in ("pident", "evalue", "score", "length", "qstart", "qend", "sstart", "send"):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        return df
    return rows

# ----------------- Summary helpers -----------------
def merge_intervals(intervals: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    if not intervals:
        return []
    iv = sorted(intervals, key=lambda x: x[0])
    merged = [iv[0]]
    for a,b in iv[1:]:
        last_a, last_b = merged[-1]
        if a <= last_b:
            merged[-1] = (last_a, max(last_b, b))
        else:
            merged.append((a,b))
    return merged

def covered_length(intervals: List[Tuple[int,int]]) -> int:
    return sum((b - a + 1) for a,b in merge_intervals(intervals))

# ----------------- Compute summary -----------------
def compute_summary_from_df(df, query_len: int, acc_len_dict: Dict[str,int]):
    """Считает сводную таблицу результатов BLAST в стиле NCBI."""
    use_pandas = pd and isinstance(df, pd.DataFrame)
    records = df.to_dict(orient="records") if use_pandas else list(df)
    groups = {}
    for r in records:
        subj = r.get("sseqid")
        if subj:
            groups.setdefault(subj, []).append(r)

    summary_rows = []
    for subj, recs in groups.items():
        scores = [Decimal(str(r.get("score") or 0.0)) for r in recs]
        max_score = max(scores) if scores else Decimal(0.0)
        total_score = sum(scores)

        evalues = [float(r.get("evalue")) for r in recs if r.get("evalue") not in (None, "")]
        if not evalues:
            min_evalue = 1.0
        else:
            min_evalue = min(evalues)
        evalue_str = "0.0" if min_evalue < 1e-180 else f"{min_evalue:.2e}"

        intervals = []
        matches_sum, aligned_sum = 0.0, 0
        for r in recs:
            qs, qe = int(r.get("qstart")), int(r.get("qend"))
            intervals.append((min(qs, qe), max(qs, qe)))
            alen = int(r.get("length") or 0)
            pid = float(r.get("pident") or 0.0)
            matches_sum += (pid / 100.0) * alen
            aligned_sum += alen

        query_cover = (covered_length(intervals) / query_len * 100.0) if query_len else 0.0
        per_ident = (matches_sum / aligned_sum * 100.0) if aligned_sum else 0.0
        acc_len = acc_len_dict.get(subj, recs[0].get("length"))

        summary_rows.append({
            "Description": subj,
            "Scientific Name": "None",
            "Max Score": int(round(max_score)),
            "Total Score": int(round(total_score)),
            "Query Cover": f"{int(round(query_cover))}%",
            "E value": evalue_str,
            "Per. Ident": f"{per_ident:.2f}%",
            "Acc. Len": acc_len
        })

    return pd.DataFrame(summary_rows) if pd else summary_rows

# ----------------- IO -----------------
def save_summary(summary, out_dir: Path, fmt: str = "both"):
    """Сохраняет сводную таблицу в CSV и/или HTML."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if pd:
        if fmt in ("csv","both"):
            summary.to_csv(out_dir / "summary.csv", index=False)
        if fmt in ("html","both"):
            summary.to_html(out_dir / "summary.html", index=False)
        print("Saved summary:", out_dir / "summary.csv", out_dir / "summary.html")
    else:
        import csv
        with (out_dir / "summary.csv").open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
            writer.writeheader()
            for r in summary:
                writer.writerow(r)

def pretty_print(summary):
    """Выводит сводную таблицу в консоль."""
    if pd:
        print(summary.to_string(index=False))
    else:
        keys = list(summary[0].keys())
        print("\t".join(keys))
        for r in summary:
            print("\t".join(str(r.get(k, "")) for k in keys))

# ----------------- Main -----------------
def main():
    parser = argparse.ArgumentParser(description="Run local blastn on query vs subjects (file or raw seq).")
    parser.add_argument("--query", required=True, help="Query FASTA path or raw sequence string")
    parser.add_argument("--subjects", required=True, nargs="+", help="One or more subject FASTA paths or raw sequence strings")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--evalue", type=float, default=1e-3)
    parser.add_argument("--max_target_seqs", type=int, default=100)
    parser.add_argument("--format", choices=["csv","html","both"], default="both", help="Save format")
    args = parser.parse_args()

    check_tools()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="blast_work_") as tmpdir:
        work = Path(tmpdir)
        query_path, _ = prepare_query_file(args.query, work)
        db_prefix, _, acc_len_dict = prepare_combined_subjects(args.subjects, work)

        out_tsv = outdir / "blast_out.tsv"
        run_blast(db_prefix, query_path, out_tsv,
                  use_db=True, threads=args.threads,
                  evalue=args.evalue, max_target_seqs=args.max_target_seqs)

        query_len = sum(len(l.strip()) for l in Path(query_path).read_text().splitlines() if not l.startswith(">"))
        results = load_results_tsv(out_tsv)
        summary = compute_summary_from_df(results, query_len, acc_len_dict)
        save_summary(summary, outdir, fmt=args.format)

        print("\n--- Summary (top results) ---")
        pretty_print(summary)

if __name__ == "__main__":
    main()
