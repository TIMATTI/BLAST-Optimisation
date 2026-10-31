import argparse
from pathlib import Path
import tempfile
from blast_utils import check_tools, run_blast
from blast_processing import prepare_query_file, prepare_combined_subjects, compute_summary_from_df
from blast_io import save_summary, pretty_print

def main():
    parser = argparse.ArgumentParser(description="Run local blastn on query vs subjects (file or raw seq).")
    parser.add_argument("--query", required=True)
    parser.add_argument("--subjects", required=True, nargs="+")
    parser.add_argument("--out", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--evalue", type=float, default=1e-3)
    parser.add_argument("--max_target_seqs", type=int, default=100)
    parser.add_argument("--format", choices=["csv","html","both"], default="both")
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
        from blast_utils import load_results_tsv
        results = load_results_tsv(out_tsv)
        summary = compute_summary_from_df(results, query_len, acc_len_dict)
        save_summary(summary, outdir, fmt=args.format)

        print("\n--- Summary ---")
        pretty_print(summary)

if __name__ == "__main__":
    main()
