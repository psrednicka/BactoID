#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
from datetime import datetime
from pathlib import Path

import xlsxwriter
from PIL import Image


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))



def read_manifest(path: Path | None) -> dict[str, str]:
    """
    Opcjonalny manifest TSV:
    Barcode<TAB>SampleName
    """
    if path is None or not path.exists():
        return {}

    mapping: dict[str, str] = {}

    with path.open(encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"Barcode", "SampleName"}

        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                "Manifest musi zawierac kolumny: Barcode i SampleName"
            )

        for row in reader:
            barcode = (row.get("Barcode") or "").strip()
            sample_name = (row.get("SampleName") or "").strip()

            if barcode and sample_name:
                mapping[barcode] = sample_name

    return mapping


def safe_float(value: str | None) -> float | None:
    if value in {"", None, "N/A"}:
        return None

    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def normalize_taxon_name(name: str) -> str:
    """
    Ujednolica nazwę do poziomu użytecznego w podsumowaniu.

    Usuwa oznaczenia szczepów, np. DSM 20284 lub ATCC 9790,
    ale zachowuje nazwę podgatunku.
    """
    cleaned = " ".join((name or "").strip().split())

    if not cleaned:
        return ""

    tokens = cleaned.split()

    if tokens[0].lower() == "candidatus" and len(tokens) >= 3:
        return " ".join(tokens[:3])

    if len(tokens) < 2:
        return cleaned

    base = tokens[:2]

    if (
        len(tokens) >= 4
        and tokens[2].lower().rstrip(".") in {"subsp", "subspecies"}
    ):
        base = tokens[:4]

    return " ".join(base)


def taxon_genus(name: str) -> str:
    tokens = name.split()

    if not tokens:
        return ""

    if tokens[0].lower() == "candidatus" and len(tokens) >= 2:
        return f"Candidatus {tokens[1]}"

    return tokens[0]


def compact_taxa(names: list[str]) -> str:
    """
    Dla jednego rodzaju:
    Lactiplantibacillus plantarum / L. pentosus

    Dla różnych rodzajów pozostawia pełne nazwy i rozdziela je znakiem +.
    """
    unique_names = list(dict.fromkeys(name for name in names if name))

    if not unique_names:
        return "No identification"

    genera = {taxon_genus(name) for name in unique_names}

    if len(genera) != 1:
        return " + ".join(unique_names)

    first = unique_names[0]
    first_tokens = first.split()
    genus = taxon_genus(first)

    compact = [first]

    for name in unique_names[1:]:
        tokens = name.split()

        if (
            genus.startswith("Candidatus ")
            or len(tokens) < 2
            or tokens[0] != first_tokens[0]
        ):
            compact.append(name)
            continue

        compact.append(f"{tokens[0][0]}. {' '.join(tokens[1:])}")

    return " / ".join(compact)


def safe_int(value: str | None) -> int:
    number = safe_float(value)

    if number is None:
        return 0

    return int(number)


def split_statuses(value: str | None) -> set[str]:
    return {
        status.strip()
        for status in (value or "").split(";")
        if status.strip()
    }


def build_cluster_records(
    sample_rows: list[dict[str, str]],
) -> list[dict]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}

    for row in sample_rows:
        key = (
            row.get("variant", ""),
            row.get("cluster_id", ""),
        )
        grouped.setdefault(key, []).append(row)

    clusters: list[dict] = []

    for (variant, cluster_id), cluster_rows in grouped.items():
        identified_rows = [
            row
            for row in cluster_rows
            if row.get("identification")
            and row.get("identification") != "No identification"
        ]

        names = sorted({
            normalize_taxon_name(
                row.get("scientific_name")
                or row.get("identification", "")
            )
            for row in identified_rows
            if normalize_taxon_name(
                row.get("scientific_name")
                or row.get("identification", "")
            )
        })

        identities = [
            value
            for value in (
                safe_float(row.get("percent_identity"))
                for row in identified_rows
            )
            if value is not None
        ]

        reference_coverages = [
            value
            for value in (
                safe_float(row.get("subject_coverage"))
                for row in identified_rows
            )
            if value is not None
        ]

        statuses: set[str] = set()

        for row in cluster_rows:
            statuses.update(split_statuses(row.get("final_status")))

        clusters.append({
            "variant": variant,
            "cluster_id": cluster_id,
            "role": cluster_rows[0].get("cluster_role", ""),
            "reads": max(
                safe_int(row.get("supporting_reads"))
                for row in cluster_rows
            ),
            "names": names,
            "genera": {taxon_genus(name) for name in names if name},
            "identities": identities,
            "reference_coverages": reference_coverages,
            "statuses": statuses,
        })

    return clusters


def select_summary_clusters(
    clusters: list[dict],
    min_reads: int = 100,
    min_fraction_of_largest: float = 0.10,
) -> list[dict]:
    """
    Do podsumowania wchodzą klastry z wiarygodną identyfikacją, które:
    - mają co najmniej min_reads,
    - stanowią co najmniej min_fraction_of_largest największego klastra.

    Jeżeli próbka ma tylko bardzo mały klaster, pokazujemy go mimo wszystko,
    ale zachowujemy ostrzeżenie o małej liczbie odczytów.
    """
    identified = [
        cluster
        for cluster in clusters
        if cluster["names"]
        and "TOO_SHORT_CONSENSUS" not in cluster["statuses"]
        and "CONCATEMER_SUSPECTED" not in cluster["statuses"]
        and "NO_IDENTIFICATION" not in cluster["statuses"]
    ]

    if not identified:
        return []

    largest = max(cluster["reads"] for cluster in identified)

    selected = [
        cluster
        for cluster in identified
        if cluster["reads"] >= min_reads
        and (
            largest == 0
            or cluster["reads"] >= largest * min_fraction_of_largest
        )
    ]

    if selected:
        return selected

    primary = [
        cluster
        for cluster in identified
        if cluster["role"] == "PRIMARY"
    ]

    if primary:
        return primary

    return [
        max(
            identified,
            key=lambda cluster: cluster["reads"],
        )
    ]


def get_primary_summary(
    rows: list[dict[str, str]],
    min_reads: int = 100,
    min_fraction_of_largest: float = 0.10,
) -> list[dict[str, str]]:
    """
    Tworzy jedno biologicznie interpretowalne podsumowanie na próbkę.

    Nie ogranicza się do klastra PRIMARY. Uwzględnia wszystkie znaczące
    klastry, dzięki czemu blisko spokrewnione gatunki mogą zostać pokazane
    jako nierozróżnialne przez 16S, zamiast arbitralnego wyboru jednego z nich.
    """
    by_sample: dict[str, list[dict[str, str]]] = {}

    for row in rows:
        by_sample.setdefault(row["sample"], []).append(row)

    summary: list[dict[str, str]] = []

    for sample, sample_rows in sorted(by_sample.items()):
        clusters = build_cluster_records(sample_rows)
        selected = select_summary_clusters(
            clusters,
            min_reads=min_reads,
            min_fraction_of_largest=min_fraction_of_largest,
        )

        if not selected:
            all_closest = sorted({
                normalize_taxon_name(row.get("closest_match", ""))
                for row in sample_rows
                if normalize_taxon_name(row.get("closest_match", ""))
            })

            summary.append({
                "sample": sample,
                "identification": "No identification",
                "reported_clusters": "0",
                "supporting_reads": "0",
                "percent_identity": "",
                "subject_coverage": "",
                "resolution": (
                    f"Closest match: {compact_taxa(all_closest)}"
                    if all_closest
                    else "No valid curated match"
                ),
                "final_status": "NO_IDENTIFICATION",
            })
            continue

        selected = sorted(
            selected,
            key=lambda cluster: (
                cluster["role"] != "PRIMARY",
                -cluster["reads"],
                cluster["variant"],
                cluster["cluster_id"],
            ),
        )

        all_names = list(dict.fromkeys(
            name
            for cluster in selected
            for name in cluster["names"]
        ))

        all_genera = {
            taxon_genus(name)
            for name in all_names
            if taxon_genus(name)
        }

        resolved_cluster_genera = {
            next(iter(cluster["genera"]))
            for cluster in selected
            if len(cluster["genera"]) == 1
        }

        identities = [
            value
            for cluster in selected
            for value in cluster["identities"]
        ]

        reference_coverages = [
            value
            for cluster in selected
            for value in cluster["reference_coverages"]
        ]

        existing_statuses: set[str] = set()

        for cluster in selected:
            existing_statuses.update(cluster["statuses"])

        ignored_summary_statuses = {
            "PASS",
            "SHORT_REFERENCE",
            "AMBIGUOUS",
        }

        warnings = sorted(
            status
            for status in existing_statuses
            if status not in ignored_summary_statuses
        )

        if len(resolved_cluster_genera) > 1:
            resolution = (
                "Substantial clusters assigned to different genera; "
                "possible mixed culture, contamination or poor isolation."
            )
            resolution_status = "POSSIBLE_MIXED_CULTURE"
        elif len(all_genera) > 1:
            resolution = (
                "One consensus has equivalent top matches from different "
                "genera; taxonomic ambiguity, not direct evidence of a mixture."
            )
            resolution_status = "TAXONOMIC_AMBIGUITY"
        elif len(all_names) > 1:
            resolution = (
                "Closely related taxa within one genus are not reliably "
                "resolved by the 16S rRNA gene. A mixture cannot be excluded."
            )
            resolution_status = "SPECIES_UNRESOLVED_BY_16S"
        else:
            resolution = "Single best 16S assignment."
            resolution_status = "PASS"

        if resolution_status != "PASS":
            warnings.append(resolution_status)

        final_status = ";".join(dict.fromkeys(warnings)) or "PASS"

        summary.append({
            "sample": sample,
            "identification": compact_taxa(all_names),
            "reported_clusters": str(len(selected)),
            "supporting_reads": str(
                sum(cluster["reads"] for cluster in selected)
            ),
            "percent_identity": (
                f"{max(identities):.3f}"
                if identities
                else ""
            ),
            "subject_coverage": (
                f"{min(reference_coverages):.3f}"
                if reference_coverages
                else ""
            ),
            "resolution": resolution,
            "final_status": final_status,
        })

    return summary


def write_value(
    worksheet,
    row_index: int,
    col_index: int,
    key: str,
    value,
    formats: dict,
) -> None:
    if key == "supporting_reads":
        if value in {"", None}:
            worksheet.write_blank(
                row_index,
                col_index,
                None,
                formats["integer"],
            )
        else:
            worksheet.write_number(
                row_index,
                col_index,
                int(float(value)),
                formats["integer"],
            )
        return

    numeric_columns = {
        "percent_identity",
        "query_coverage",
        "query_length",
        "subject_length",
        "alignment_length",
        "subject_coverage",
        "closest_match_identity",
        "closest_match_query_coverage",
        "closest_match_subject_coverage",
    }

    if key in numeric_columns:
        number = safe_float(value)
        cell_format = (
            formats["integer"]
            if key in {
                "query_length",
                "subject_length",
                "alignment_length",
            }
            else formats["decimal"]
        )

        if number is None:
            worksheet.write_blank(
                row_index,
                col_index,
                None,
                cell_format,
            )
        else:
            worksheet.write_number(
                row_index,
                col_index,
                number,
                cell_format,
            )
        return

    worksheet.write(
        row_index,
        col_index,
        value or "",
        formats["text"],
    )


def add_status_formatting(
    worksheet,
    first_row: int,
    last_row: int,
    status_column: int,
    formats: dict,
) -> None:
    if last_row < first_row:
        return

    # Najpoważniejszy status ma pierwszeństwo.
    conditions = [
        ("NO_IDENTIFICATION", formats["error"], True),
        ("TOO_SHORT_CONSENSUS", formats["error"], True),
        ("CONCATEMER_SUSPECTED", formats["error"], True),
        ("VERY_LOW_READ_COUNT", formats["warning"], False),
        ("LOW_READ_COUNT", formats["warning"], False),
        ("LOW_REFERENCE_COVERAGE", formats["warning"], False),
        ("LOW_QUERY_COVERAGE", formats["warning"], False),
        ("POSSIBLE_MIXED_CULTURE", formats["warning"], False),
        ("SPECIES_UNRESOLVED_BY_16S", formats["warning"], False),
        ("TAXONOMIC_AMBIGUITY", formats["warning"], False),
        ("AMBIGUOUS", formats["warning"], False),
        ("PASS", formats["pass"], False),
    ]

    for value, cell_format, stop_if_true in conditions:
        worksheet.conditional_format(
            first_row,
            status_column,
            last_row,
            status_column,
            {
                "type": "text",
                "criteria": "containing",
                "value": value,
                "format": cell_format,
                "stop_if_true": stop_if_true,
            },
        )


def calculate_logo_placement(
    logo_path: Path,
    panel_width_px: int = 650,
    panel_height_px: int = 270,
    margin_px: int = 8,
    fill_ratio: float = 0.96,
) -> dict[str, float | int]:
    """
    Dopasuj logo do panelu G1:J7 bez deformowania proporcji.

    Logo może być zarówno pomniejszane, jak i powiększane.
    Poprzednia wersja ograniczała skalę do maksymalnie 1.0,
    przez co małe obrazy nie były powiększane.
    """
    with Image.open(logo_path) as image:
        width_px, height_px = image.size

    if width_px <= 0 or height_px <= 0:
        raise ValueError(f"Nieprawidłowe wymiary logo: {logo_path}")

    available_width = (
        panel_width_px - (2 * margin_px)
    ) * fill_ratio
    available_height = (
        panel_height_px - (2 * margin_px)
    ) * fill_ratio

    # Bez ograniczenia do 1.0: logo może zostać powiększone.
    scale = min(
        available_width / width_px,
        available_height / height_px,
    )

    rendered_width = width_px * scale
    rendered_height = height_px * scale

    return {
        "x_scale": scale,
        "y_scale": scale,
        "x_offset": max(
            margin_px,
            int((panel_width_px - rendered_width) / 2),
        ),
        "y_offset": max(
            margin_px,
            int((panel_height_px - rendered_height) / 2),
        ),
        "object_position": 1,
    }


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--manifest",
        default="",
        help="Optional TSV file with columns Barcode and SampleName",
    )
    parser.add_argument("--pipeline-version", required=True)
    parser.add_argument("--ngspeciesid-version", required=True)
    parser.add_argument("--blast-version", required=True)
    parser.add_argument("--database-name", required=True)
    parser.add_argument("--database-path", required=True)
    parser.add_argument("--logo", default="")
    parser.add_argument("--min-read-length", type=int, default=1000)
    parser.add_argument("--max-read-length", type=int, default=2000)
    parser.add_argument("--summary-min-reads", type=int, default=100)
    parser.add_argument(
        "--summary-min-cluster-fraction",
        type=float,
        default=0.10,
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    logo_path = Path(args.logo) if args.logo else None
    database_path = Path(args.database_path)
    manifest_path = Path(args.manifest) if args.manifest else None

    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = read_tsv(input_path)

    sample_mapping = read_manifest(manifest_path)

    for row in rows:
        original_sample = row.get("sample", "")
        row["sample"] = sample_mapping.get(
            original_sample,
            original_sample,
        )

    if not rows:
        raise RuntimeError(
            f"Pusty plik wejściowy: {input_path}"
        )

    summary_rows = get_primary_summary(
        rows,
        min_reads=args.summary_min_reads,
        min_fraction_of_largest=args.summary_min_cluster_fraction,
    )

    database_files = sorted(
        database_path.parent.glob(database_path.name + ".*")
    )

    if database_files:
        database_date = max(
            path.stat().st_mtime
            for path in database_files
        )
        database_date_text = datetime.fromtimestamp(
            database_date
        ).strftime("%Y-%m-%d %H:%M")
    else:
        database_date_text = "Unknown"

    generated = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    workbook = xlsxwriter.Workbook(str(output_path))

    workbook.set_properties({
        "title": "BactoID genetic identification report",
        "subject": "Identification of bacterial isolates",
        "author": "Laboratory of Molecular Biology",
        "company": "IBPRS-PIB",
        "comments": f"BactoID version {args.pipeline_version}",
    })

    formats = {
        "title": workbook.add_format({
            "bold": True,
            "font_size": 19,
            "font_color": "#FFFFFF",
            "bg_color": "#123B5D",
            "align": "left",
            "valign": "vcenter",
        }),
        "subtitle": workbook.add_format({
            "bold": True,
            "font_size": 11,
            "font_color": "#D9EAF2",
            "bg_color": "#123B5D",
            "align": "left",
            "valign": "vcenter",
        }),
        "workflow_label": workbook.add_format({
            "bold": True,
            "font_color": "#A9DDE0",
            "bg_color": "#123B5D",
            "align": "left",
            "valign": "vcenter",
        }),
        "workflow_value": workbook.add_format({
            "font_color": "#FFFFFF",
            "bg_color": "#123B5D",
            "align": "left",
            "valign": "vcenter",
            "text_wrap": True,
        }),
        "logo_panel": workbook.add_format({
            "bg_color": "#FFFFFF",
            "border": 1,
            "border_color": "#D7E1E8",
            "align": "center",
            "valign": "vcenter",
        }),
        "section": workbook.add_format({
            "bold": True,
            "font_size": 12,
            "font_color": "#FFFFFF",
            "bg_color": "#007C83",
            "border": 1,
            "align": "left",
            "valign": "vcenter",
        }),
        "header": workbook.add_format({
            "bold": True,
            "font_color": "#FFFFFF",
            "bg_color": "#007C83",
            "border": 1,
            "align": "center",
            "valign": "vcenter",
            "text_wrap": True,
        }),
        "text": workbook.add_format({
            "border": 1,
            "valign": "top",
            "text_wrap": True,
        }),
        "integer": workbook.add_format({
            "border": 1,
            "align": "center",
            "valign": "vcenter",
            "num_format": "0",
        }),
        "decimal": workbook.add_format({
            "border": 1,
            "align": "center",
            "valign": "vcenter",
            "num_format": "0.000",
        }),
        "label": workbook.add_format({
            "bold": True,
            "bg_color": "#DDEBF7",
            "border": 1,
        }),
        "value": workbook.add_format({
            "border": 1,
            "text_wrap": True,
        }),
        "pass": workbook.add_format({
            "bg_color": "#C6EFCE",
            "font_color": "#006100",
        }),
        "warning": workbook.add_format({
            "bg_color": "#FFF2CC",
            "font_color": "#7F6000",
        }),
        "error": workbook.add_format({
            "bg_color": "#F4CCCC",
            "font_color": "#9C0006",
        }),
        "note": workbook.add_format({
            "font_color": "#666666",
            "italic": True,
            "text_wrap": True,
        }),
    }

    # ========================================================
    # SUMMARY
    # ========================================================

    summary_ws = workbook.add_worksheet("Summary")
    summary_ws.hide_gridlines(2)
    summary_ws.set_tab_color("#007C83")
    summary_ws.set_zoom(85)

    # --------------------------------------------------------
    # Nagłówek:
    # - po lewej A1:F7: informacje o workflow,
    # - po prawej G1:J7: duże logo BactoID.
    # --------------------------------------------------------

    summary_ws.set_column("A:B", 17)
    summary_ws.set_column("C:F", 20)
    summary_ws.set_column("G:J", 22)

    # Wyższy nagłówek daje logo realnie więcej miejsca.
    for row_number in range(0, 7):
        summary_ws.set_row(row_number, 28)

    # Lewy panel – workflow i wersjonowanie.
    summary_ws.merge_range(
        "A1:F1",
        "BactoID workflow report",
        formats["title"],
    )
    summary_ws.merge_range(
        "A2:F2",
        "Genetic identification of bacterial isolates from ONT data",
        formats["subtitle"],
    )

    workflow_information = [
        ("Generated", generated),
        (
            "Software",
            (
                f"BactoID {args.pipeline_version}  |  "
                f"NGSpeciesID {args.ngspeciesid_version}  |  "
                f"BLAST+ {args.blast_version}"
            ),
        ),
        ("Database", args.database_name),
        ("Database files", database_date_text),
        ("Database path", str(database_path)),
    ]

    for row_index, (label, value) in enumerate(
        workflow_information,
        start=2,
    ):
        summary_ws.merge_range(
            row_index,
            0,
            row_index,
            1,
            label,
            formats["workflow_label"],
        )
        summary_ws.merge_range(
            row_index,
            2,
            row_index,
            5,
            value,
            formats["workflow_value"],
        )

    # Prawy panel – wyłącznie logo.
    if logo_path and logo_path.exists():
        summary_ws.merge_range(
            "G1:J7",
            "",
            formats["logo_panel"],
        )
        placement = calculate_logo_placement(logo_path)
        summary_ws.insert_image(
            "G1",
            str(logo_path),
            placement,
        )
    else:
        summary_ws.merge_range(
            "G1:J7",
            "BactoID",
            formats["logo_panel"],
        )

    # Odstęp między nagłówkiem a tabelą.
    summary_ws.set_row(7, 8)

    summary_ws.merge_range(
        "A9:H9",
        "Identification summary",
        formats["section"],
    )

    summary_columns = [
        ("sample", "Sample"),
        ("identification", "Identification based on 16S"),
        ("reported_clusters", "Reported clusters"),
        ("supporting_reads", "Reads in reported clusters"),
        ("percent_identity", "Best identity [%]"),
        ("subject_coverage", "Minimum reference coverage [%]"),
        ("resolution", "Interpretation"),
        ("final_status", "Status"),
    ]

    # Wiersz 10 w Excelu = indeks 9.
    header_row = 9

    for column_index, (_, label) in enumerate(summary_columns):
        summary_ws.write(
            header_row,
            column_index,
            label,
            formats["header"],
        )

    for row_index, row in enumerate(
        summary_rows,
        start=header_row + 1,
    ):
        for column_index, (key, _) in enumerate(summary_columns):
            write_value(
                summary_ws,
                row_index,
                column_index,
                key,
                row.get(key, ""),
                formats,
            )

    summary_last_row = header_row + len(summary_rows)

    summary_ws.autofilter(
        header_row,
        0,
        summary_last_row,
        len(summary_columns) - 1,
    )

    summary_ws.freeze_panes(header_row + 1, 0)

    # Szerokości tabeli.
    summary_ws.set_column("A:A", 15)
    summary_ws.set_column("B:B", 42)
    summary_ws.set_column("C:C", 17)
    summary_ws.set_column("D:D", 23)
    summary_ws.set_column("E:F", 22)
    summary_ws.set_column("G:G", 62)
    summary_ws.set_column("H:H", 38)

    add_status_formatting(
        summary_ws,
        header_row + 1,
        summary_last_row,
        7,
        formats,
    )

    # ========================================================
    # ALL CLUSTERS
    # ========================================================

    all_ws = workbook.add_worksheet("All clusters")
    all_ws.hide_gridlines(2)
    all_ws.freeze_panes(1, 0)
    all_ws.set_tab_color("#4F81BD")

    detail_columns = [
        ("sample", "Sample"),
        ("variant", "Variant"),
        ("cluster_id", "Cluster ID"),
        ("cluster_role", "Cluster role"),
        ("supporting_reads", "Supporting reads"),
        ("read_status", "Read status"),
        ("identification", "Identification"),
        ("scientific_name", "Scientific name"),
        ("percent_identity", "Identity [%]"),
        ("query_coverage", "Query coverage [%]"),
        ("query_length", "Query length [bp]"),
        ("subject_length", "Reference length [bp]"),
        ("alignment_length", "Alignment length [bp]"),
        ("subject_coverage", "Reference coverage [%]"),
        ("accession", "Accession"),
        ("closest_match", "Closest match"),
        (
            "closest_match_identity",
            "Closest identity [%]",
        ),
        (
            "closest_match_query_coverage",
            "Closest query coverage [%]",
        ),
        (
            "closest_match_subject_coverage",
            "Closest reference coverage [%]",
        ),
        (
            "closest_match_accession",
            "Closest accession",
        ),
        ("amplicon_length_status", "Amplicon length status"),
        ("coverage_status", "Coverage status"),
        ("ambiguity_status", "Ambiguity"),
        ("final_status", "Final status"),
    ]

    for column_index, (_, label) in enumerate(detail_columns):
        all_ws.write(
            0,
            column_index,
            label,
            formats["header"],
        )

    for row_index, row in enumerate(rows, start=1):
        for column_index, (key, _) in enumerate(detail_columns):
            write_value(
                all_ws,
                row_index,
                column_index,
                key,
                row.get(key, ""),
                formats,
            )

    all_ws.autofilter(
        0,
        0,
        len(rows),
        len(detail_columns) - 1,
    )

    all_ws.set_column("A:A", 14)
    all_ws.set_column("B:E", 14)
    all_ws.set_column("F:F", 24)
    all_ws.set_column("G:H", 29)
    all_ws.set_column("I:J", 17)
    all_ws.set_column("K:M", 18)
    all_ws.set_column("N:N", 21)
    all_ws.set_column("O:O", 16)
    all_ws.set_column("P:P", 29)
    all_ws.set_column("Q:S", 20)
    all_ws.set_column("T:T", 18)
    all_ws.set_column("U:W", 24)
    all_ws.set_column("X:X", 48)

    add_status_formatting(
        all_ws,
        1,
        len(rows),
        23,
        formats,
    )

    # ========================================================
    # INTERPRETATION
    # ========================================================

    interpretation_ws = workbook.add_worksheet("Interpretation")
    interpretation_ws.hide_gridlines(2)
    interpretation_ws.set_tab_color("#F4B183")

    interpretation_ws.merge_range(
        "A1:D1",
        "BactoID interpretation rules",
        formats["title"],
    )

    interpretation_ws.set_column("A:A", 27)
    interpretation_ws.set_column("B:B", 20)
    interpretation_ws.set_column("C:C", 25)
    interpretation_ws.set_column("D:D", 75)

    interpretation_columns = [
        "Criterion",
        "Threshold",
        "Status",
        "Interpretation",
    ]

    for index, label in enumerate(interpretation_columns):
        interpretation_ws.write(
            2,
            index,
            label,
            formats["header"],
        )

    interpretation_rows = [
        [
            "Consensus length",
            f"{args.min_read_length}–{args.max_read_length} bp",
            "PASS",
            (
                "Expected interval for a bacterial full-length 16S amplicon "
                "in BactoID. Reads outside this interval are removed before "
                "NGSpeciesID consensus generation."
            ),
        ],
        [
            "Consensus length",
            f"<{args.min_read_length} bp",
            "TOO_SHORT_CONSENSUS",
            (
                "The sequence is too short for full-length 16S "
                "identification and is not accepted as a final result."
            ),
        ],
        [
            "Consensus length",
            f">{args.max_read_length} bp",
            "CONCATEMER_SUSPECTED",
            (
                "The sequence is longer than the expected 16S amplicon. "
                "A tandem amplicon, concatemer or chimeric sequence is "
                "suspected and no final identification is reported."
            ),
        ],
        [
            "Percent identity",
            "≥97.0%",
            "IDENTIFIED",
            (
                "The best curated BLAST match satisfies the minimum "
                "identity threshold and the consensus length is valid."
            ),
        ],
        [
            "Percent identity",
            "<97.0%",
            "NO_IDENTIFICATION",
            (
                "No taxonomic identification is reported. The closest "
                "curated match is retained for technical review."
            ),
        ],
        [
            "Reference coverage",
            "≥95%",
            "PASS / SHORT_REFERENCE",
            (
                "Most of the reference sequence is covered by the alignment. "
                "Reference coverage is shown in the Summary sheet."
            ),
        ],
        [
            "Reference coverage",
            "<95%",
            "LOW_REFERENCE_COVERAGE",
            (
                "A substantial part of the reference is not covered and the "
                "result requires cautious interpretation."
            ),
        ],
        [
            "Supporting reads",
            "≥100",
            "PASS",
            "Adequate read support for consensus generation.",
        ],
        [
            "Supporting reads",
            "20–99",
            "LOW_READ_COUNT",
            "Identification is reported with a low-read support warning.",
        ],
        [
            "Supporting reads",
            "1–19",
            "VERY_LOW_READ_COUNT",
            (
                "Identification is based on very limited read support and "
                "requires confirmation."
            ),
        ],
        [
            "Several reportable taxa",
            "Same genus",
            "SPECIES_UNRESOLVED_BY_16S",
            (
                "All substantial clusters belong to one genus, but more than "
                "one species or subspecies is supported. Report all plausible "
                "taxa separated by a slash. This is not, by itself, proof of "
                "a mixed culture; 16S may lack discriminatory power."
            ),
        ],
        [
            "Several reportable clusters",
            "Different genera",
            "POSSIBLE_MIXED_CULTURE",
            (
                "Substantial clusters assigned independently to different "
                "genera suggest a mixed culture, contamination, barcode "
                "cross-talk or inadequate colony isolation."
            ),
        ],
        [
            "Equivalent top hits",
            "Different genera in one consensus",
            "TAXONOMIC_AMBIGUITY",
            (
                "One consensus has equivalent top matches from different "
                "genera. This is taxonomic ambiguity and is not direct "
                "evidence that the sample contains several organisms."
            ),
        ],
        [
            "Equivalent top hits",
            "Several species",
            "AMBIGUOUS",
            (
                "Several curated references share the same maximum percent "
                "identity. Full-length 16S may not provide unambiguous "
                "species resolution."
            ),
        ],
    ]

    for row_index, row in enumerate(
        interpretation_rows,
        start=3,
    ):
        for column_index, value in enumerate(row):
            interpretation_ws.write(
                row_index,
                column_index,
                value,
                formats["text"],
            )

    # ========================================================
    # ABOUT
    # ========================================================

    about_ws = workbook.add_worksheet("About")
    about_ws.hide_gridlines(2)
    about_ws.set_tab_color("#A5A5A5")

    about_ws.merge_range(
        "A1:D1",
        "BactoID run information",
        formats["title"],
    )

    about_ws.set_column("A:A", 27)
    about_ws.set_column("B:B", 55)
    about_ws.set_column("C:D", 20)

    run_information = [
        ("Generated", generated),
        ("BactoID version", args.pipeline_version),
        ("NGSpeciesID version", args.ngspeciesid_version),
        ("BLAST version", args.blast_version),
        ("Database", args.database_name),
        ("Database path", str(database_path)),
        ("Database files modified", database_date_text),
        ("Number of samples", str(len(summary_rows))),
        ("Number of reported rows", str(len(rows))),
        (
            "Minimum identification threshold",
            "97.0% identity",
        ),
        (
            "Accepted read length",
            f"{args.min_read_length}–{args.max_read_length} bp",
        ),
        (
            "Minimum recommended reference coverage",
            "95%",
        ),
        (
            "Length filtering",
            (
                "Reads shorter than the minimum and longer than the maximum "
                "are removed before NGSpeciesID."
            ),
        ),
    ]

    for row_index, (label, value) in enumerate(
        run_information,
        start=2,
    ):
        about_ws.write(
            row_index,
            0,
            label,
            formats["label"],
        )
        about_ws.write(
            row_index,
            1,
            value,
            formats["value"],
        )

    about_ws.write(
        len(run_information) + 4,
        0,
        (
            "BactoID performs consensus-based identification "
            "of bacterial isolates using full-length 16S rRNA "
            "Oxford Nanopore sequencing reads."
        ),
        formats["note"],
    )

    workbook.close()


if __name__ == "__main__":
    main()
