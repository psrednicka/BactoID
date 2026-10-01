#!/usr/bin/env python3

import argparse
import csv
import re
import shutil
from pathlib import Path


HEADER_PATTERN = re.compile(
    r"consensus_cl_id_(?P<cluster>\d+)_total_supporting_reads_(?P<reads>\d+)"
)


def read_fasta_header(path: Path) -> tuple[str, int]:
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().strip()

    if not header.startswith(">"):
        raise ValueError(f"Nieprawidłowy FASTA: {path}")

    match = HEADER_PATTERN.search(header)

    if not match:
        raise ValueError(
            f"Nie można odczytać klastra i liczby odczytów z nagłówka: {header}"
        )

    return match.group("cluster"), int(match.group("reads"))


def read_status(reads: int) -> str:
    if reads >= 100:
        return "PASS"
    if reads >= 20:
        return "LOW_READ_COUNT"
    return "VERY_LOW_READ_COUNT"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    summary_path = Path(args.summary)

    consensus_files = sorted(
        input_dir.glob("racon_cl_id_*/consensus.fasta")
    )

    if not consensus_files:
        raise RuntimeError(
            f"Nie znaleziono konsensusów po Medace w: {input_dir}"
        )

    clusters = []

    for fasta in consensus_files:
        cluster_id, reads = read_fasta_header(fasta)

        clusters.append({
            "cluster_id": cluster_id,
            "supporting_reads": reads,
            "source_fasta": fasta,
        })

    clusters.sort(
        key=lambda item: item["supporting_reads"],
        reverse=True,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []

    for variant, cluster in enumerate(clusters, start=1):
        destination = output_dir / f"cluster_{cluster['cluster_id']}.fasta"
        shutil.copy2(cluster["source_fasta"], destination)

        rows.append({
            "variant": variant,
            "cluster_id": cluster["cluster_id"],
            "cluster_role": "PRIMARY" if variant == 1 else "SECONDARY",
            "supporting_reads": cluster["supporting_reads"],
            "read_status": read_status(cluster["supporting_reads"]),
            "source_fasta": str(cluster["source_fasta"]),
            "blast_input_fasta": str(destination),
        })

    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "variant",
                "cluster_id",
                "cluster_role",
                "supporting_reads",
                "read_status",
                "source_fasta",
                "blast_input_fasta",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
