#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
from pathlib import Path
from typing import Any


BLAST_COLUMNS = [
    "qseqid",
    "saccver",
    "pident",
    "qcovs",
    "qlen",
    "slen",
    "length",
    "qstart",
    "qend",
    "sstart",
    "send",
    "bitscore",
    "sscinames",
    "stitle",
]

OUTPUT_COLUMNS = [
    "sample",
    "variant",
    "cluster_id",
    "cluster_role",
    "supporting_reads",
    "read_status",
    "identification",
    "scientific_name",
    "percent_identity",
    "query_coverage",
    "query_length",
    "subject_length",
    "alignment_length",
    "subject_coverage",
    "accession",
    "closest_match",
    "closest_match_identity",
    "closest_match_query_coverage",
    "closest_match_subject_coverage",
    "closest_match_accession",
    "amplicon_length_status",
    "coverage_status",
    "ambiguity_status",
    "final_status",
    "title",
]


def safe_float(value: Any) -> float | None:
    if value in {None, "", "N/A"}:
        return None

    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_blast(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []

    hits: list[dict[str, str]] = []

    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue

            values = line.rstrip("\n").split("\t")

            if len(values) != len(BLAST_COLUMNS):
                raise RuntimeError(
                    f"{path}: wiersz {line_number} ma {len(values)} kolumn, "
                    f"a oczekiwano {len(BLAST_COLUMNS)}. Usuń stare wyniki "
                    "BLAST i uruchom pipeline ponownie po zmianie outfmt."
                )

            hits.append(dict(zip(BLAST_COLUMNS, values)))

    return hits


def get_cluster_value(row: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = row.get(key, "")
        if value not in {None, ""}:
            return str(value)
    return ""


def find_top_hits_file(
    blast_dir: Path,
    variant: str,
    cluster_id: str,
) -> Path:
    candidates: list[Path] = []

    for value in [variant, cluster_id]:
        if not value:
            continue

        if value.startswith("cluster_"):
            candidates.append(blast_dir / f"{value}_top.tsv")
        else:
            candidates.append(blast_dir / f"cluster_{value}_top.tsv")
            candidates.append(blast_dir / f"{value}_top.tsv")

    for candidate in candidates:
        if candidate.exists():
            return candidate

    unique_candidates = list(dict.fromkeys(candidates))

    if unique_candidates:
        return unique_candidates[0]

    raise RuntimeError(
        f"Nie można ustalić nazwy pliku BLAST dla klastra: "
        f"{variant=} {cluster_id=}"
    )


def calculate_subject_coverage(hit: dict[str, str]) -> float | None:
    subject_length = safe_float(hit.get("slen"))
    subject_start = safe_float(hit.get("sstart"))
    subject_end = safe_float(hit.get("send"))

    if (
        subject_length is None
        or subject_length <= 0
        or subject_start is None
        or subject_end is None
    ):
        return None

    subject_span = abs(subject_end - subject_start) + 1
    return min(100.0, (subject_span / subject_length) * 100.0)


def classify_amplicon_length(
    hit: dict[str, str],
    min_query_length: int,
    max_query_length: int,
) -> str:
    query_length = safe_float(hit.get("qlen"))

    if query_length is None:
        return "LENGTH_NOT_DETERMINED"

    if query_length < min_query_length:
        return "TOO_SHORT_CONSENSUS"

    if query_length > max_query_length:
        return "CONCATEMER_SUSPECTED"

    return "PASS"


def classify_coverage(
    hit: dict[str, str],
    min_query_coverage: float,
    min_subject_coverage: float,
) -> tuple[str, float | None]:
    query_coverage = safe_float(hit.get("qcovs"))
    query_length = safe_float(hit.get("qlen"))
    subject_length = safe_float(hit.get("slen"))
    subject_coverage = calculate_subject_coverage(hit)

    if query_coverage is None or subject_coverage is None:
        return "NOT_DETERMINED", subject_coverage

    if (
        query_coverage >= min_query_coverage
        and subject_coverage >= min_subject_coverage
    ):
        return "PASS", subject_coverage

    if (
        query_length is not None
        and subject_length is not None
        and subject_length < query_length
        and subject_coverage >= min_subject_coverage
    ):
        return "SHORT_REFERENCE", subject_coverage

    if subject_coverage < min_subject_coverage:
        return "LOW_REFERENCE_COVERAGE", subject_coverage

    return "LOW_QUERY_COVERAGE", subject_coverage


def build_final_status(
    identified: bool,
    read_status: str,
    amplicon_length_status: str,
    coverage_status: str,
    ambiguity_status: str,
) -> str:
    statuses: list[str] = []

    if not identified:
        statuses.append("NO_IDENTIFICATION")

    if read_status and read_status != "PASS":
        statuses.append(read_status)

    if amplicon_length_status != "PASS":
        statuses.append(amplicon_length_status)

    if coverage_status in {
        "LOW_QUERY_COVERAGE",
        "LOW_REFERENCE_COVERAGE",
        "NOT_DETERMINED",
    }:
        statuses.append(coverage_status)

    if ambiguity_status == "AMBIGUOUS":
        statuses.append(ambiguity_status)

    if not statuses:
        if coverage_status == "SHORT_REFERENCE":
            return "PASS;SHORT_REFERENCE"
        return "PASS"

    if coverage_status == "SHORT_REFERENCE":
        statuses.append("SHORT_REFERENCE")

    return ";".join(dict.fromkeys(statuses))


def make_empty_result(
    sample: str,
    cluster: dict[str, str],
) -> dict[str, str]:
    variant = get_cluster_value(cluster, "variant")
    cluster_id = get_cluster_value(cluster, "cluster_id", "cluster")
    cluster_role = get_cluster_value(cluster, "cluster_role", "role")
    supporting_reads = get_cluster_value(cluster, "supporting_reads", "reads")
    read_status = get_cluster_value(cluster, "read_status")

    result = {column: "" for column in OUTPUT_COLUMNS}
    result.update({
        "sample": sample,
        "variant": variant,
        "cluster_id": cluster_id,
        "cluster_role": cluster_role,
        "supporting_reads": supporting_reads,
        "read_status": read_status,
        "identification": "No identification",
        "amplicon_length_status": "NOT_DETERMINED",
        "coverage_status": "NO_BLAST_RESULT",
        "ambiguity_status": "NOT_APPLICABLE",
        "final_status": "NO_IDENTIFICATION",
    })
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--blast-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-identity", type=float, default=97.0)
    parser.add_argument("--min-query-coverage", type=float, default=90.0)
    parser.add_argument("--min-subject-coverage", type=float, default=95.0)
    parser.add_argument("--min-query-length", type=int, default=1000)
    parser.add_argument("--max-query-length", type=int, default=2000)
    args = parser.parse_args()

    clusters_path = Path(args.clusters)
    blast_dir = Path(args.blast_dir)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    clusters = read_tsv(clusters_path)
    sample = clusters_path.parent.name
    output_rows: list[dict[str, str]] = []

    for cluster in clusters:
        variant = get_cluster_value(cluster, "variant")
        cluster_id = get_cluster_value(cluster, "cluster_id", "cluster")
        cluster_role = get_cluster_value(cluster, "cluster_role", "role")
        supporting_reads = get_cluster_value(
            cluster,
            "supporting_reads",
            "reads",
        )
        read_status = get_cluster_value(cluster, "read_status")

        top_hits_path = find_top_hits_file(
            blast_dir,
            variant,
            cluster_id,
        )
        hits = read_blast(top_hits_path)

        if not hits:
            output_rows.append(make_empty_result(sample, cluster))
            continue

        scientific_names = sorted({
            hit.get("sscinames", "").strip()
            for hit in hits
            if hit.get("sscinames", "").strip()
        })

        ambiguity_status = (
            "AMBIGUOUS"
            if len(scientific_names) > 1
            else "UNAMBIGUOUS"
        )

        for hit in hits:
            identity = safe_float(hit.get("pident"))
            amplicon_length_status = classify_amplicon_length(
                hit,
                min_query_length=args.min_query_length,
                max_query_length=args.max_query_length,
            )

            # Nie raportujemy formalnej identyfikacji, jeśli konsensus nie ma
            # biologicznie wiarygodnej długości pełnego 16S.
            identified = (
                identity is not None
                and identity >= args.min_identity
                and amplicon_length_status == "PASS"
            )

            coverage_status, subject_coverage = classify_coverage(
                hit,
                min_query_coverage=args.min_query_coverage,
                min_subject_coverage=args.min_subject_coverage,
            )

            if amplicon_length_status != "PASS":
                coverage_status = "NOT_EVALUATED_INVALID_LENGTH"

            scientific_name = hit.get("sscinames", "").strip()
            subject_coverage_text = (
                f"{subject_coverage:.3f}"
                if subject_coverage is not None
                else ""
            )

            row = {column: "" for column in OUTPUT_COLUMNS}
            row.update({
                "sample": sample,
                "variant": variant,
                "cluster_id": cluster_id,
                "cluster_role": cluster_role,
                "supporting_reads": supporting_reads,
                "read_status": read_status,
                "query_length": hit.get("qlen", ""),
                "subject_length": hit.get("slen", ""),
                "alignment_length": hit.get("length", ""),
                "subject_coverage": subject_coverage_text,
                "amplicon_length_status": amplicon_length_status,
                "coverage_status": coverage_status,
                "ambiguity_status": ambiguity_status,
                "final_status": build_final_status(
                    identified=identified,
                    read_status=read_status,
                    amplicon_length_status=amplicon_length_status,
                    coverage_status=coverage_status,
                    ambiguity_status=ambiguity_status,
                ),
                "title": hit.get("stitle", ""),
            })

            if identified:
                row.update({
                    "identification": scientific_name,
                    "scientific_name": scientific_name,
                    "percent_identity": hit.get("pident", ""),
                    "query_coverage": hit.get("qcovs", ""),
                    "accession": hit.get("saccver", ""),
                })
            else:
                row.update({
                    "identification": "No identification",
                    "closest_match": scientific_name,
                    "closest_match_identity": hit.get("pident", ""),
                    "closest_match_query_coverage": hit.get("qcovs", ""),
                    "closest_match_subject_coverage": subject_coverage_text,
                    "closest_match_accession": hit.get("saccver", ""),
                })

            output_rows.append(row)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
