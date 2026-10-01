#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
from pathlib import Path


def iter_fastq(path: Path):
    with path.open(encoding="utf-8") as handle:
        record_number = 0

        while True:
            header = handle.readline()
            if not header:
                break

            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            record_number += 1

            if not sequence or not plus or not quality:
                raise RuntimeError(
                    f"Niepełny rekord FASTQ nr {record_number} w pliku {path}"
                )

            if not header.startswith("@"):
                raise RuntimeError(
                    f"Nieprawidłowy nagłówek FASTQ w rekordzie {record_number}: "
                    f"{header.rstrip()}"
                )

            if not plus.startswith("+"):
                raise RuntimeError(
                    f"Brak linii '+' w rekordzie FASTQ nr {record_number}"
                )

            sequence_text = sequence.rstrip("\r\n")
            quality_text = quality.rstrip("\r\n")

            if len(sequence_text) != len(quality_text):
                raise RuntimeError(
                    f"Długość sekwencji i jakości różni się w rekordzie "
                    f"FASTQ nr {record_number}"
                )

            yield (
                header,
                sequence,
                plus,
                quality,
                len(sequence_text),
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Filtruje odczyty FASTQ według długości dla pełnej długości "
            "ampliconu bakteryjnego 16S."
        )
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--min-length", type=int, default=1000)
    parser.add_argument("--max-length", type=int, default=2000)
    args = parser.parse_args()

    if args.min_length <= 0:
        raise ValueError("--min-length musi być większe od zera")

    if args.max_length < args.min_length:
        raise ValueError("--max-length nie może być mniejsze niż --min-length")

    input_path = Path(args.input)
    output_path = Path(args.output)
    stats_path = Path(args.stats)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    stats_path.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    kept = 0
    too_short = 0
    too_long = 0
    min_observed = None
    max_observed = None

    with output_path.open("w", encoding="utf-8", newline="") as output_handle:
        for header, sequence, plus, quality, length in iter_fastq(input_path):
            total += 1
            min_observed = length if min_observed is None else min(min_observed, length)
            max_observed = length if max_observed is None else max(max_observed, length)

            if length < args.min_length:
                too_short += 1
                continue

            if length > args.max_length:
                too_long += 1
                continue

            kept += 1
            output_handle.write(header)
            output_handle.write(sequence)
            output_handle.write(plus)
            output_handle.write(quality)

    if total == 0:
        raise RuntimeError(f"Brak odczytów w pliku wejściowym: {input_path}")

    if kept == 0:
        raise RuntimeError(
            f"Po filtracji {args.min_length}-{args.max_length} bp nie pozostał "
            f"żaden odczyt: {input_path}"
        )

    with stats_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "input",
            "min_length",
            "max_length",
            "total_reads",
            "kept_reads",
            "too_short",
            "too_long",
            "kept_percent",
            "min_observed_length",
            "max_observed_length",
        ])
        writer.writerow([
            str(input_path),
            args.min_length,
            args.max_length,
            total,
            kept,
            too_short,
            too_long,
            f"{(kept / total) * 100:.3f}",
            min_observed,
            max_observed,
        ])


if __name__ == "__main__":
    main()
