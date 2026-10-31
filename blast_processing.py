# blast_processing.py
from pathlib import Path
from typing import List, Tuple, Dict
from decimal import Decimal
import math

# безопасный импорт pandas — модуль может работать и без pandas
try:
    import pandas as pd
except Exception:
    pd = None

# импортируем утилиты (убедись, что blast_utils.py находится рядом и экспортирует эти функции)
from blast_utils import fasta_wrap, write_fasta_from_seq, makeblastdb, load_results_tsv

# ----------------- FASTA парсер -----------------
def parse_fasta_entries(text: str) -> List[Tuple[str, str]]:
    entries = []
    cur_h = None
    cur_seq_parts = []
    for line in text.splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if cur_h is not None:
                entries.append((cur_h, "".join(cur_seq_parts)))
            cur_h = line[1:].strip()
            cur_seq_parts = []
        else:
            cur_seq_parts.append(line.strip())
    if cur_h is not None:
        entries.append((cur_h, "".join(cur_seq_parts)))
    return entries

# ----------------- Prepare query and subjects -----------------
def prepare_query_file(query_arg: str, workdir: Path) -> Tuple[str, int]:
    p = Path(query_arg)
    if p.exists():
        text = p.read_text()
        if text.lstrip().startswith(">"):
            entries = parse_fasta_entries(text)
            header, seq = entries[0]
        else:
            seq = "".join(text.split())
        return str(p), len(seq)
    qf = workdir / "query.fasta"
    write_fasta_from_seq(query_arg, qf, header="query")
    return str(qf), len("".join(query_arg.split()))

def prepare_combined_subjects(subject_args: List[str], workdir: Path) -> Tuple[str, List[Path], Dict[str,int]]:
    created = []
    seqs = []
    acc_len_dict: Dict[str,int] = {}

    for i, subj in enumerate(subject_args):
        p = Path(subj)
        if p.exists():
            text = p.read_text()
            if text.lstrip().startswith(">"):
                entries = parse_fasta_entries(text)
                for (header, seq) in entries:
                    header_id = header.split()[0]
                    seqs.append(f">{header_id}\n{fasta_wrap(seq)}")
                    acc_len_dict[header_id] = len(seq)
            else:
                header_id = f"subject_file_{i}"
                seq = "".join(text.split())
                seqs.append(f">{header_id}\n{fasta_wrap(seq)}")
                acc_len_dict[header_id] = len(seq)
        else:
            header_id = f"subject_{i}"
            seq = subj
            seqs.append(f">{header_id}\n{fasta_wrap(seq)}")
            acc_len_dict[header_id] = len(seq)

    combined_fasta = workdir / "combined_subjects.fasta"
    with combined_fasta.open("w", encoding="utf-8") as fh:
        for s in seqs:
            fh.write(s.strip() + "\n")
    created.append(combined_fasta)

    db_prefix = workdir / "subjects_db"
    makeblastdb(combined_fasta, db_prefix)
    created.append(db_prefix)
    return str(db_prefix), created, acc_len_dict

# ----------------- Summary helpers -----------------
def merge_intervals(intervals: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    if not intervals:
        return []
    iv = sorted(intervals, key=lambda x: x[0])
    merged = [iv[0]]
    for a, b in iv[1:]:
        last_a, last_b = merged[-1]
        if a <= last_b:
            merged[-1] = (last_a, max(last_b, b))
        else:
            merged.append((a, b))
    return merged

def covered_length(intervals: List[Tuple[int,int]]) -> int:
    return sum((b - a + 1) for a, b in merge_intervals(intervals))

# ----------------- Compute summary -----------------
def compute_summary_from_df(df, query_len: int, acc_len_dict: Dict[str,int]):
    """
    Robust summary computation:
     - prefer bitscore then score,
     - prefer slen, then acc_len_dict, then length as fallback,
     - filter NaN/inf everywhere,
     - return pandas DataFrame if available, else list of dicts.
    """
    use_pandas = (pd is not None) and isinstance(df, pd.DataFrame)
    records = df.to_dict(orient="records") if use_pandas else list(df)

    # collect by subject id
    groups: Dict[str, List[dict]] = {}
    for r in records:
        subj = r.get("sseqid")
        if subj:
            groups.setdefault(subj, []).append(r)

    # small diagnostic: if fields missing, we note it (non-verbose)
    if not groups:
        print("DEBUG: no groups parsed from BLAST output (check OUTFMT alignment).", file=sys.stderr)

    summary_rows = []
    for subj, recs in groups.items():
        # collect scores (prefer bitscore -> score), filter NaN/inf
        scores: List[Decimal] = []
        for r in recs:
            sval_num = None
            for fld in ("bitscore", "score"):
                candidate = r.get(fld)
                if candidate not in (None, ""):
                    try:
                        v = float(candidate)
                        if math.isfinite(v):
                            sval_num = v
                            break
                    except Exception:
                        continue
            if sval_num is not None:
                try:
                    scores.append(Decimal(str(sval_num)))
                except Exception:
                    scores.append(Decimal(sval_num))

        max_score = max(scores) if scores else Decimal(0)
        total_score = sum(scores, Decimal(0))

        # evalue
        evalues: List[float] = []
        for r in recs:
            ev_raw = r.get("evalue")
            if ev_raw not in (None, ""):
                try:
                    ev = float(ev_raw)
                    if math.isfinite(ev):
                        evalues.append(ev)
                except Exception:
                    continue
        min_evalue = min(evalues) if evalues else 1.0
        evalue_str = "0.0" if min_evalue < 1e-180 else f"{min_evalue:.2e}"

        # coverage & percent identity
        intervals: List[Tuple[int,int]] = []
        matches_sum = 0.0
        aligned_sum = 0
        for r in recs:
            try:
                qs = int(float(r.get("qstart")))
                qe = int(float(r.get("qend")))
                intervals.append((min(qs, qe), max(qs, qe)))
            except Exception:
                pass
            try:
                alen = int(float(r.get("length") or 0))
            except Exception:
                alen = 0
            try:
                pid = float(r.get("pident") or 0.0)
            except Exception:
                pid = 0.0
            matches_sum += (pid / 100.0) * alen
            aligned_sum += alen

        query_cover = (covered_length(intervals) / query_len * 100.0) if query_len and intervals else 0.0
        per_ident = (matches_sum / aligned_sum * 100.0) if aligned_sum else 0.0

        # Acc. Len: prefer slen -> acc_len_dict -> fallback length of first alignment
        acc_len = None
        for r in recs:
            sval = r.get("slen")
            if sval not in (None, ""):
                try:
                    svalf = float(sval)
                    if math.isfinite(svalf):
                        acc_len = int(svalf)
                        break
                except Exception:
                    pass
        if acc_len is None:
            acc_len = acc_len_dict.get(subj)
        if acc_len is None:
            try:
                acc_len = int(float(recs[0].get("length") or 0))
            except Exception:
                acc_len = 0

        # safe convert Decimal->int
        try:
            max_score_int = int(round(float(max_score)))
        except Exception:
            max_score_int = 0
        try:
            total_score_int = int(round(float(total_score)))
        except Exception:
            total_score_int = 0

        summary_rows.append({
            "Description": subj,
            "Scientific Name": "None",
            "Max Score": max_score_int,
            "Total Score": total_score_int,
            "Query Cover": f"{int(round(query_cover))}%",
            "E value": evalue_str,
            "Per. Ident": f"{per_ident:.2f}%",
            "Acc. Len": acc_len
        })

    return pd.DataFrame(summary_rows) if (pd is not None) else summary_rows
