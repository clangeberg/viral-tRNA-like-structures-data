#!/usr/bin/env python3
"""Recalculate alignment entropy and insertion-frequency tables.

The calculations used for the structure-mapped panels are intentionally kept
dependency-free and use the checked-in Stockholm alignments and residue maps.
Henikoff position-based sequence weights are calculated over all alignment
columns. Shannon entropy is calculated from the weighted, non-gap A/C/G/U
frequencies and normalized by log(4), so 0 is invariant and 1 is maximally
diverse. Insertion frequencies follow the reference-gap-boundary definition
described in the manuscript.
"""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean


ROOT = Path(__file__).resolve().parents[1]
ALIGNMENTS = ROOT / "alignments"
DATA = ROOT / "data"
BASES = "ACGU"
NO_INSERTION_COLOR = "#D3D3D3"
MAX_INSERTION_COLOR = "#B7FF00"

PANEL_CONFIGS = (
    {
        "class": "TLSVal",
        "alignment": "TLSVal.sto",
        "mapping": "TLSVal_structure_mapping.csv",
        "entropy": "TLSVal_alignment_entropy.csv",
    },
    {
        "class": "TLSHis",
        "alignment": "TLSHis.sto",
        "mapping": "TLSHis_structure_mapping.csv",
        "entropy": "TLSHis_alignment_entropy.csv",
    },
    {
        "class": "TLSTyr Bromovirus (E/B3 subclass)",
        "alignment": "TLSTyr_Bromovirus.sto",
        "mapping": "TLSTyr_Bromovirus_structure_mapping.csv",
        "entropy": "TLSTyr_Bromovirus_alignment_entropy.csv",
    },
)


def parse_stockholm(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = defaultdict(list)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped == "//":
                continue
            fields = stripped.split()
            if len(fields) >= 2:
                records[fields[0]].append(fields[1])
    joined = {
        name: "".join(parts).upper().replace("T", "U")
        for name, parts in records.items()
    }
    lengths = {len(sequence) for sequence in joined.values()}
    if not joined or len(lengths) != 1:
        raise ValueError(f"Inconsistent alignment lengths in {path}: {lengths}")
    return joined


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def henikoff_metrics(records: dict[str, str]) -> tuple[list[dict], dict[str, float]]:
    names = list(records)
    sequences = [records[name] for name in names]
    width = len(sequences[0])
    raw_weights = [0.0] * len(names)

    for column in range(width):
        column_bases = [
            sequence[column] if sequence[column] in BASES else None
            for sequence in sequences
        ]
        counts = Counter(base for base in column_bases if base is not None)
        diversity = len(counts)
        if not diversity:
            continue
        for index, base in enumerate(column_bases):
            if base is not None:
                raw_weights[index] += 1.0 / (diversity * counts[base])

    total_weight = sum(raw_weights)
    weights = (
        [weight / total_weight for weight in raw_weights]
        if total_weight
        else [1.0 / len(names)] * len(names)
    )
    named_weights = dict(zip(names, weights))

    rows: list[dict] = []
    for column in range(width):
        weighted_counts = {base: 0.0 for base in BASES}
        raw_counts: Counter[str] = Counter()
        for index, sequence in enumerate(sequences):
            base = sequence[column]
            if base in BASES:
                weighted_counts[base] += weights[index]
                raw_counts[base] += 1
            else:
                raw_counts["gap_or_other"] += 1
        occupancy = sum(weighted_counts.values())
        if occupancy:
            frequencies = {
                base: weighted_counts[base] / occupancy for base in BASES
            }
            entropy = -sum(
                frequency * math.log(frequency)
                for frequency in frequencies.values()
                if frequency > 0
            ) / math.log(4.0)
            consensus = max(BASES, key=lambda base: frequencies[base])
            conservation = 1.0 - entropy
        else:
            frequencies = {base: 0.0 for base in BASES}
            entropy = math.nan
            conservation = math.nan
            consensus = "-"
        rows.append(
            {
                "sto_column_1based": column + 1,
                "consensus_base": consensus,
                "normalized_shannon_entropy": entropy,
                "conservation_score_1_minus_entropy": conservation,
                "henikoff_weighted_occupancy": occupancy,
                **{f"weighted_frequency_{base}": frequencies[base] for base in BASES},
                **{f"raw_count_{base}": raw_counts[base] for base in BASES},
                "raw_count_gap_or_other": raw_counts["gap_or_other"],
            }
        )
    return rows, named_weights


def reference_column_maps(aligned_reference: str) -> tuple[dict[int, int], dict[int, int]]:
    position_to_column: dict[int, int] = {}
    column_to_position: dict[int, int] = {}
    reference_position = 0
    for column, character in enumerate(aligned_reference):
        if character in BASES:
            reference_position += 1
            position_to_column[reference_position] = column
            column_to_position[column] = reference_position
    return position_to_column, column_to_position


def reference_gap_blocks(aligned_reference: str):
    start = None
    for column, character in enumerate(aligned_reference):
        if character not in BASES and start is None:
            start = column
        elif character in BASES and start is not None:
            yield start, column
            start = None
    if start is not None:
        yield start, len(aligned_reference)


def rgb(hex_color: str) -> tuple[int, int, int]:
    return tuple(int(hex_color[index:index + 2], 16) for index in (1, 3, 5))


def insertion_color(prevalence: float, scale_maximum: float) -> str:
    fraction = 0.0 if scale_maximum <= 0 else min(1.0, prevalence / scale_maximum)
    low = rgb(NO_INSERTION_COLOR)
    high = rgb(MAX_INSERTION_COLOR)
    values = [round(a + (b - a) * fraction) for a, b in zip(low, high)]
    return "#" + "".join(f"{value:02X}" for value in values)


def calculate_insertions(
    class_name: str,
    alignment_file: str,
    records: dict[str, str],
    weights: dict[str, float],
    mapping_rows: list[dict[str, str]],
) -> tuple[list[dict], dict[int, float], dict[int, list[str]]]:
    reference_name = mapping_rows[0]["reference_record"]
    reference = records[reference_name]
    _, column_to_position = reference_column_maps(reference)
    mapping_by_reference_position = {
        int(row["reference_residue_1based"]): row for row in mapping_rows
    }
    boundaries: list[dict] = []
    score_by_residue: dict[int, float] = {}
    boundaries_by_residue: dict[int, list[str]] = defaultdict(list)

    for start, end in reference_gap_blocks(reference):
        left_column = start - 1
        right_column = end
        if left_column not in column_to_position or right_column not in column_to_position:
            continue
        left_position = column_to_position[left_column]
        right_position = column_to_position[right_column]
        evaluable: list[str] = []
        inserted: list[str] = []
        insertion_lengths: list[int] = []
        inserted_all = 0
        for name, sequence in records.items():
            insertion_length = sum(base in BASES for base in sequence[start:end])
            if insertion_length:
                inserted_all += 1
            if sequence[left_column] not in BASES or sequence[right_column] not in BASES:
                continue
            evaluable.append(name)
            if insertion_length:
                inserted.append(name)
                insertion_lengths.append(insertion_length)
        if not evaluable or not inserted:
            continue

        boundary_id = f"{alignment_file.removesuffix('.sto')}_I{len(boundaries) + 1:02d}"
        raw_prevalence = len(inserted) / len(evaluable)
        evaluable_weight = sum(weights[name] for name in evaluable)
        weighted_prevalence = (
            sum(weights[name] for name in inserted) / evaluable_weight
            if evaluable_weight else 0.0
        )
        left_mapping = mapping_by_reference_position.get(left_position)
        right_mapping = mapping_by_reference_position.get(right_position)
        colored_residues: list[int] = []
        for flank in (left_mapping, right_mapping):
            if not flank or flank["resolved_in_model_1"].lower() != "true":
                continue
            residue = int(flank["pdb_label_seq_id"])
            colored_residues.append(residue)
            score_by_residue[residue] = max(score_by_residue.get(residue, 0.0), raw_prevalence)
            boundaries_by_residue[residue].append(boundary_id)

        boundaries.append(
            {
                "class_or_subclass": class_name,
                "boundary_id": boundary_id,
                "alignment_file": alignment_file,
                "reference_record": reference_name,
                "sto_block_start_1based": start + 1,
                "sto_block_end_1based": end,
                "aligned_block_width": end - start,
                "left_reference_residue_1based": left_position,
                "right_reference_residue_1based": right_position,
                "left_reference_base": reference[left_column],
                "right_reference_base": reference[right_column],
                "alignment_record_count": len(records),
                "evaluable_record_count": len(evaluable),
                "inserted_record_count_evaluable": len(inserted),
                "inserted_record_count_all": inserted_all,
                "raw_insertion_prevalence_evaluable": raw_prevalence,
                "raw_insertion_prevalence_all": inserted_all / len(records),
                "henikoff_weighted_insertion_prevalence_evaluable": weighted_prevalence,
                "mean_insertion_length_nt": mean(insertion_lengths),
                "maximum_insertion_length_nt": max(insertion_lengths),
                "left_pdb_label_seq_id": left_mapping["pdb_label_seq_id"] if left_mapping else "NA",
                "right_pdb_label_seq_id": right_mapping["pdb_label_seq_id"] if right_mapping else "NA",
                "left_flank_resolved": left_mapping["resolved_in_model_1"] if left_mapping else "False",
                "right_flank_resolved": right_mapping["resolved_in_model_1"] if right_mapping else "False",
                "colored_pdb_label_seq_ids": ";".join(map(str, sorted(set(colored_residues)))) or "NA",
                "inserted_records": ";".join(inserted),
            }
        )
    return boundaries, score_by_residue, dict(boundaries_by_residue)


def main() -> None:
    all_alignment_files = sorted(ALIGNMENTS.glob("*.sto"))
    count_rows = []
    for path in all_alignment_files:
        records = parse_stockholm(path)
        count_rows.append(
            {
                "alignment": path.name,
                "sequence_records": len(records),
                "alignment_columns": len(next(iter(records.values()))),
            }
        )
    write_csv(
        DATA / "sequence_counts.csv",
        count_rows,
        ["alignment", "sequence_records", "alignment_columns"],
    )

    all_boundaries: list[dict] = []
    panel_results: list[dict] = []
    for config in PANEL_CONFIGS:
        records = parse_stockholm(ALIGNMENTS / config["alignment"])
        entropy_rows, weights = henikoff_metrics(records)
        entropy_fields = list(entropy_rows[0])
        write_csv(DATA / "entropy" / config["entropy"], entropy_rows, entropy_fields)

        mapping_path = DATA / "mappings" / config["mapping"]
        mapping_rows = read_csv(mapping_path)
        reference_name = mapping_rows[0]["reference_record"]
        position_to_column, _ = reference_column_maps(records[reference_name])
        structure_entropy_rows: list[dict] = []
        for mapping in mapping_rows:
            reference_position = int(mapping["reference_residue_1based"])
            column = position_to_column[reference_position]
            structure_entropy_rows.append(
                {
                    **mapping,
                    **entropy_rows[column],
                }
            )
        structure_entropy_filename = config["mapping"].replace(
            "_structure_mapping.csv", "_structure_entropy_mapping.csv"
        )
        write_csv(
            DATA / "mappings" / structure_entropy_filename,
            structure_entropy_rows,
            list(structure_entropy_rows[0]),
        )

        boundaries, score_by_residue, boundary_ids = calculate_insertions(
            config["class"], config["alignment"], records, weights, mapping_rows
        )
        all_boundaries.extend(boundaries)
        panel_results.append(
            {
                "class": config["class"],
                "mapping_rows": mapping_rows,
                "score_by_residue": score_by_residue,
                "boundary_ids": boundary_ids,
            }
        )

    scale_maximum = max(
        float(row["raw_insertion_prevalence_evaluable"])
        for row in all_boundaries
    )
    for row in all_boundaries:
        prevalence = float(row["raw_insertion_prevalence_evaluable"])
        row["normalized_global_color_intensity"] = prevalence / scale_maximum
        row["boundary_color_hex"] = insertion_color(prevalence, scale_maximum)
    write_csv(
        DATA / "insertions" / "TLS_insertion_boundaries.csv",
        all_boundaries,
        list(all_boundaries[0]),
    )

    residue_rows: list[dict] = []
    for result in panel_results:
        mapping = result["mapping_rows"][0]
        for residue, prevalence in sorted(result["score_by_residue"].items()):
            residue_rows.append(
                {
                    "class_or_subclass": result["class"],
                    "pdb_id": mapping["pdb_id"],
                    "pdb_label_chain": mapping["pdb_label_chain"],
                    "pdb_label_seq_id": residue,
                    "maximum_flanking_insertion_prevalence": prevalence,
                    "normalized_global_color_intensity": prevalence / scale_maximum,
                    "color_hex": insertion_color(prevalence, scale_maximum),
                    "boundary_ids": ";".join(result["boundary_ids"][residue]),
                }
            )
    write_csv(
        DATA / "insertions" / "TLS_insertion_residue_colors.csv",
        residue_rows,
        list(residue_rows[0]),
    )

    scale_rows = []
    for tick in range(11):
        intensity = tick / 10
        prevalence = scale_maximum * intensity
        scale_rows.append(
            {
                "scale_id": "InsertionFrequencyLime_GlobalLinear",
                "tick_index": tick,
                "normalized_intensity": intensity,
                "insertion_prevalence": prevalence,
                "insertion_prevalence_percent": prevalence * 100,
                "color_hex": insertion_color(prevalence, scale_maximum),
                "scale_minimum_prevalence": 0.0,
                "scale_maximum_prevalence": scale_maximum,
                "scale_minimum_color_hex": NO_INSERTION_COLOR,
                "scale_maximum_color_hex": MAX_INSERTION_COLOR,
                "normalization": "linear; shared globally across the three structure-mapped panels",
                "score_definition": "records with at least one nucleotide in a reference-gap block divided by evaluable records retaining both reference flanks",
            }
        )
    write_csv(
        DATA / "insertions" / "TLS_insertion_color_scale_metadata.csv",
        scale_rows,
        list(scale_rows[0]),
    )

    print(f"Recalculated {len(PANEL_CONFIGS)} entropy panels")
    print(f"Global insertion-frequency maximum: {scale_maximum:.9f} ({scale_maximum:.1%})")
    for config in PANEL_CONFIGS:
        rows = [row for row in all_boundaries if row["class_or_subclass"] == config["class"]]
        print(
            f"{config['class']}: {len(rows)} insertion boundaries; "
            f"maximum {max(float(row['raw_insertion_prevalence_evaluable']) for row in rows):.1%}"
        )


if __name__ == "__main__":
    main()
