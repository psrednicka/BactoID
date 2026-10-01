#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    classification_dir = Path(args.classification_dir)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    files = sorted(classification_dir.glob("barcode*.tsv"))

    if not files:
        raise RuntimeError(
            f"Brak plików klasyfikacji w: {classification_dir}"
        )

    rows = []
    fieldnames = None

    for path in files:
        sample = path.stem

        with path.open(encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")

            if fieldnames is None:
                fieldnames = ["sample"] + reader.fieldnames

            for row in reader:
                rows.append({
                    "sample": sample,
                    **row,
                })

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
