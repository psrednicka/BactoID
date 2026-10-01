#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path


def parse_identity(line: str) -> float:
    columns = line.rstrip("\n").split("\t")

    if len(columns) < 3:
        raise ValueError(
            "Nieprawidłowy wiersz BLAST: oczekiwano co najmniej 3 kolumn."
        )

    try:
        return float(columns[2])
    except ValueError as exc:
        raise ValueError(
            f"Nieprawidłowa wartość pident w wierszu: {line.rstrip()}"
        ) from exc


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Zachowuje wszystkie hity BLAST z najwyższym percent identity. "
            "Próg identyfikacji jest oceniany później w classify_hits.py."
        )
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--min-identity",
        type=float,
        default=97.0,
        help=(
            "Zachowany dla zgodności z pipeline; nie usuwa trafień poniżej "
            "progu, ponieważ są potrzebne jako closest match."
        ),
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not input_path.exists() or input_path.stat().st_size == 0:
        output_path.write_text("", encoding="utf-8")
        return

    lines = [
        line
        for line in input_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]

    if not lines:
        output_path.write_text("", encoding="utf-8")
        return

    identities = [parse_identity(line) for line in lines]
    maximum_identity = max(identities)

    selected = [
        line
        for line, identity in zip(lines, identities)
        if abs(identity - maximum_identity) < 1e-9
    ]

    output_path.write_text(
        "\n".join(selected) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
