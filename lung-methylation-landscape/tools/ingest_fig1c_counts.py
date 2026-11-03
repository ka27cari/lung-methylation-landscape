#!/usr/bin/env python3
"""Append the local Figure 1c aggregate to a figure-labelled Source Data directory.

Unknown screening eligibility is not drawn in the original Figure 1c code and
is therefore excluded. This tool never writes individual-level study records.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import csv
import hashlib
from pathlib import Path

COLUMNS = ["eligibility", "diagnosis_interval", "ldct_category", "n"]
DATASET = "ldct_screening_timing_counts"


def ingest(input_csv: Path, source_dir: Path) -> tuple[int, int]:
    with input_csv.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != COLUMNS:
            raise ValueError("Unexpected columns in aggregate file")
        rows = list(reader)
    keys = [(r["eligibility"], r["diagnosis_interval"], r["ldct_category"]) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate aggregate groups")
    totals = defaultdict(int)
    for row in rows:
        if row["eligibility"] not in {"Eligible", "Ineligible", "Unknown"}:
            raise ValueError("Unexpected eligibility value")
        if row["diagnosis_interval"] not in {"<3 yrs", ">=3 yrs"}:
            raise ValueError("Unexpected timing category")
        if row["ldct_category"] not in {"Normal", "Nodule present"}:
            raise ValueError("Unexpected LDCT category")
        try:
            n = int(row["n"])
        except ValueError as exc:
            raise ValueError("Non-integer count") from exc
        if n < 0:
            raise ValueError("Negative count")
        totals[row["eligibility"]] += n
    if totals["Eligible"] != 9 or totals["Ineligible"] != 45:
        raise ValueError("Figure 1c author-reported 9/45 groups do not match")
    displayed = [r for r in rows if r["eligibility"] != "Unknown"]
    displayed.sort(key=lambda r: (
        {"Ineligible": 0, "Eligible": 1}[r["eligibility"]],
        {">=3 yrs": 0, "<3 yrs": 1}[r["diagnosis_interval"]],
        {"Normal": 0, "Nodule present": 1}[r["ldct_category"]],
    ))
    target = source_dir / "Fig1" / f"{DATASET}.csv"
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(displayed)
    manifest = source_dir / "DATASETS.csv"
    if manifest.is_file():
        with manifest.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            columns = reader.fieldnames
            index = list(reader)
        needed = ["figure", "panel", "dataset", "file", "rows", "upstream", "sha256"]
        if columns != needed:
            raise ValueError("Unexpected Source Data manifest schema")
        index = [r for r in index if not (r["figure"] == "Fig1" and r["dataset"] == DATASET)]
        index.append({"figure": "Fig1", "panel": "c", "dataset": DATASET,
                      "file": f"Fig1/{DATASET}.csv", "rows": str(len(displayed)),
                      "upstream": "author local LDCT report aggregate; Figure 1c displayed groups",
                      "sha256": hashlib.sha256(target.read_bytes()).hexdigest()})
        fig1 = sorted((r for r in index if r["figure"] == "Fig1"), key=lambda r: r["panel"])
        others = [r for r in index if r["figure"] != "Fig1"]
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=needed, lineterminator="\n")
            writer.writeheader()
            writer.writerows(fig1 + others)
    return len(displayed), sum(int(r["n"]) for r in displayed)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("figure_source_data_dir", type=Path)
    args = parser.parse_args()
    rows, people = ingest(args.input_csv, args.figure_source_data_dir)
    print(f"Figure 1c: {rows} aggregate cells, {people} displayed cases")


if __name__ == "__main__":
    main()
