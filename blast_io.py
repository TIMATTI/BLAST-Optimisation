from pathlib import Path
from typing import List, Dict

def save_summary(summary, out_dir: Path, fmt: str = "both"):
    import pandas as pd
    out_dir.mkdir(parents=True, exist_ok=True)
    if pd and hasattr(summary, "to_csv"):
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
    import pandas as pd
    if pd and hasattr(summary, "to_string"):
        print(summary.to_string(index=False))
    else:
        keys = list(summary[0].keys())
        print("\t".join(keys))
        for r in summary:
            print("\t".join(str(r.get(k, "")) for k in keys))
