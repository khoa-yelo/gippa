"""Run selection proposal from ``selection_rules.tsv`` and the per-dataset ``files.tsv``.

Stdlib only:  python3 -I validation/campaign/select_runs.py validation/campaign

Every rule's ``pattern`` is a regular expression matched in full against the raw file names of
``<accession>/files.tsv``; ``NA`` writes one placeholder row (a decision with no file yet). A rule
that matches no file is an error (exit 1), so a typo cannot silently drop runs. Writes
``run_selection.tsv`` and prints the selected GB per dataset and priority.
"""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

COLUMNS = ("accession", "arm", "priority", "file_name", "bytes", "checksum", "raw_format",
           "acquisition", "label", "status", "reason")  # fmt: skip


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main(root: Path) -> int:
    out, errors = [], []
    gb: dict[tuple[str, str], float] = defaultdict(float)
    for rule in read(root / "selection_rules.tsv"):
        acc = rule["accession"]
        base = {k: rule[k] for k in ("accession", "arm", "priority", "label", "status", "reason")}
        if rule["pattern"] == "NA":
            out.append({**base, **dict.fromkeys(("file_name", "bytes", "checksum", "raw_format",
                                                 "acquisition"), "NA")})  # fmt: skip
            continue
        pat = re.compile(rule["pattern"])
        hits = [f for f in read(root / acc / "files.tsv")
                if f["file_type"] == "raw" and pat.fullmatch(f["file_name"])]  # fmt: skip
        if not hits:
            errors.append(f"{acc} {rule['arm']}: no file matches {rule['pattern']}")
        for f in hits:
            out.append({**base, "file_name": f["file_name"], "bytes": f["bytes"],
                        "checksum": f["checksum"], "raw_format": f["raw_format"],
                        "acquisition": f["acquisition"]})  # fmt: skip
            if f["bytes"].isdigit():
                gb[(acc, rule["priority"])] += int(f["bytes"]) / 1e9
    with (root / "run_selection.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    for (acc, prio), v in sorted(gb.items()):
        print(f"{acc}\tpriority {prio}\t{v:.1f} GB")
    for e in errors:
        print(e, file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1] if len(sys.argv) > 1 else "validation/campaign")))
