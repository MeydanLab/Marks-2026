#!/usr/bin/env python3
"""Collapse adjacent peak-comparison CSV rows after threshold filtering.

Rules implemented:
- Dataset columns are all columns after 'Nucleotide Window'.
- User can use all datasets, include only named datasets, or exclude named datasets.
- A row is kept if at least N selected datasets have value > cutoff.
- Within each gene, adjacent qualifying codon rows are collapsed into one region.
  Once a run of consecutive peak codons starts, it keeps extending until the next
  codon is not a qualifying peak row. For collapsed groups:
  - keep metadata from the representative regional peak codon
  - sum dataset values across grouped rows

How this script runs, at a high level:
1) Read a peak-comparison CSV that already contains per-codon dataset values.
2) Decide which dataset columns are in play.
3) Keep only rows where enough selected datasets exceed the cutoff.
4) Within each gene, merge neighboring codon rows that behave like one region.
5) Write a simplified CSV where each merged region is represented by one row.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

GENE_COL = "Gene"
PEAK_REGION_COL = "Peak Region"
CODON_COL = "Codon Number"
NUC_WINDOW_COL = "Nucleotide Window"


def parse_args() -> argparse.Namespace:
    # This script is meant to be used from the command line after a peak CSV
    # already exists. The arguments control:
    # - which CSV to read
    # - where to write the collapsed result
    # - which numeric threshold defines a peak
    # - how many datasets must support a row
    # - whether to include all dataset columns or only a chosen subset
    parser = argparse.ArgumentParser(description="Filter and collapse peak comparison CSV rows.")
    parser.add_argument("input_csv", help="Input CSV path")
    parser.add_argument(
        "-o",
        "--output-csv",
        default=None,
        help="Output CSV path (default: <input_stem>_collapsed.csv)",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        required=True,
        help="Threshold a dataset value must exceed to count as observed",
    )
    parser.add_argument(
        "--min-datasets",
        type=int,
        required=True,
        help="Minimum selected datasets that must exceed cutoff for a row to qualify",
    )

    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--include-datasets",
        default=None,
        help="Comma-separated dataset names to include (all others ignored)",
    )
    selection.add_argument(
        "--exclude-datasets",
        default=None,
        help="Comma-separated dataset names to exclude (use all others)",
    )

    return parser.parse_args()


def split_csv_list(value: str | None) -> List[str]:
    # Include/exclude dataset names are supplied as comma-separated text on the
    # command line. This helper converts that text into a clean Python list.
    if not value:
        return []
    return [x.strip() for x in value.split(",") if x.strip()]


def parse_float(value: str) -> float:
    # Peak-comparison CSVs may contain blank cells. Here blanks are treated as
    # zero so missing values do not break later threshold comparisons.
    value = (value or "").strip()
    if value == "":
        return 0.0
    try:
        return float(value)
    except ValueError:
        raise ValueError(f"Non-numeric dataset value encountered: {value!r}")


def parse_codon(value: str) -> int:
    # Codon numbers should be integers, but CSVs sometimes store them in a form
    # like "25.0". Parsing through float first makes the script tolerant of that.
    text = (value or "").strip()
    try:
        return int(float(text))
    except ValueError:
        raise ValueError(f"Invalid codon number: {value!r}")


def dataset_columns(headers: Sequence[str]) -> Tuple[int, List[str]]:
    # The script assumes all dataset signal columns come after the
    # "Nucleotide Window" column. That lets it work on peak-comparison CSVs with
    # varying numbers of datasets without hardcoding their names.
    if NUC_WINDOW_COL not in headers:
        raise ValueError(f"Required column not found: {NUC_WINDOW_COL!r}")
    nuc_idx = headers.index(NUC_WINDOW_COL)
    cols = list(headers[nuc_idx + 1 :])
    if not cols:
        raise ValueError("No dataset columns found after 'Nucleotide Window'.")
    return nuc_idx, cols


def select_datasets(all_datasets: Sequence[str], include: List[str], exclude: List[str]) -> List[str]:
    # This block resolves which dataset columns will participate in filtering
    # and collapsing:
    # - include list: only those named datasets
    # - exclude list: all datasets except those named
    # - neither: use every dataset column
    all_set = set(all_datasets)

    if include:
        missing = [d for d in include if d not in all_set]
        if missing:
            raise ValueError(f"Included dataset(s) not found in input CSV: {', '.join(missing)}")
        return include

    if exclude:
        missing = [d for d in exclude if d not in all_set]
        if missing:
            raise ValueError(f"Excluded dataset(s) not found in input CSV: {', '.join(missing)}")
        return [d for d in all_datasets if d not in set(exclude)]

    return list(all_datasets)


def qualifying_rows(
    rows: Iterable[Dict[str, str]],
    all_datasets: Sequence[str],
    selected_datasets: Sequence[str],
    cutoff: float,
    min_datasets: int,
) -> List[Dict[str, object]]:
    # This is the first major transformation step.
    #
    # For each input row:
    # 1) read the gene/codon metadata
    # 2) parse every dataset value into numbers
    # 3) record which selected datasets exceed the cutoff
    # 4) keep the row only if at least min_datasets support it
    #
    # The output keeps both the original row text and parsed numeric metadata so
    # later steps can collapse regions without reparsing the CSV.
    kept: List[Dict[str, object]] = []

    for raw in rows:
        if GENE_COL not in raw:
            raise ValueError(f"Required column not found: {GENE_COL!r}")
        if CODON_COL not in raw:
            raise ValueError(f"Required column not found: {CODON_COL!r}")

        gene = (raw.get(GENE_COL) or "").strip()
        codon_num = parse_codon(raw.get(CODON_COL, ""))

        # Convert all dataset columns to floats once so every downstream step can
        # work with numeric values instead of repeatedly parsing strings.
        dataset_vals: Dict[str, float] = {}
        for ds in all_datasets:
            dataset_vals[ds] = parse_float(raw.get(ds, ""))

        # Track exactly which selected datasets are above threshold for this
        # codon. This is still useful for reporting/inspection, even though the
        # collapsing step now groups any consecutive qualifying codons together.
        above: Set[str] = set()
        for ds in selected_datasets:
            value = dataset_vals[ds]
            if value > cutoff:
                above.add(ds)

        if len(above) >= min_datasets:
            kept.append(
                {
                    "raw": dict(raw),
                    "gene": gene,
                    "codon": codon_num,
                    "dataset_vals": dataset_vals,
                    "above": above,
                }
            )

    return kept


def select_representative_row(
    group: List[Dict[str, object]], selected_datasets: Sequence[str]
) -> Dict[str, object]:
    # After several codons are merged into one region, the script still needs
    # one row to donate the displayed metadata such as codon number and local
    # sequence window.
    #
    # The representative is chosen by a simple voting scheme:
    # - each selected dataset "votes" for the codon where that dataset is
    #   strongest within the group
    # - the codon with the most votes wins
    # - ties are broken by strongest signal, then by earlier codon position
    if len(group) == 1:
        return group[0]

    vote_counts: Dict[int, int] = {}
    vote_strength: Dict[int, float] = {}

    for ds in selected_datasets:
        winner = max(group, key=lambda row: float(row["dataset_vals"][ds]))
        codon = int(winner["codon"])
        value = float(winner["dataset_vals"][ds])
        vote_counts[codon] = vote_counts.get(codon, 0) + 1
        previous_strength = vote_strength.get(codon, float("-inf"))
        if value > previous_strength:
            vote_strength[codon] = value

    max_votes = max(vote_counts.values())
    candidate_codons = [codon for codon, count in vote_counts.items() if count == max_votes]

    if len(candidate_codons) == 1:
        chosen_codon = candidate_codons[0]
    else:
        chosen_codon = max(candidate_codons, key=lambda codon: (vote_strength[codon], -codon))

    return next(row for row in group if int(row["codon"]) == chosen_codon)


def collapse_gene_rows(
    rows_for_gene: List[Dict[str, object]],
    all_datasets: Sequence[str],
    selected_datasets: Sequence[str],
) -> List[Dict[str, object]]:
    # This is the second major transformation step.
    #
    # It operates one gene at a time and merges rows whenever the codons are
    # consecutive after threshold filtering.
    #
    # In other words: once one codon qualifies as a peak, any immediately
    # neighboring qualifying codon is treated as part of the same regional peak,
    # regardless of which replicate or dataset contributed the signal.
    if not rows_for_gene:
        return []

    # Sorting by codon number guarantees that "adjacent" means neighboring codons
    # in the transcript, not just neighboring rows in the original file.
    rows_sorted = sorted(rows_for_gene, key=lambda r: int(r["codon"]))
    groups: List[List[Dict[str, object]]] = []

    # Build groups of rows that belong to the same regional peak.
    for row in rows_sorted:
        if not groups:
            groups.append([row])
            continue

        prev_row = groups[-1][-1]
        current_codon = int(row["codon"])
        prev_codon = int(prev_row["codon"])

        adjacent = current_codon == prev_codon + 1

        if adjacent:
            groups[-1].append(row)
        else:
            groups.append([row])

    # Once groups are defined, turn each group into one output row by:
    # - choosing a representative codon for metadata
    # - summing dataset values across the whole region
    collapsed: List[Dict[str, object]] = []
    for group in groups:
        representative = select_representative_row(group, selected_datasets)
        out_raw = dict(representative["raw"])
        start_codon = int(group[0]["codon"])
        end_codon = int(group[-1]["codon"])

        sums = {ds: 0.0 for ds in all_datasets}
        for row in group:
            for ds in all_datasets:
                sums[ds] += float(row["dataset_vals"][ds])

        out_raw[PEAK_REGION_COL] = f"{start_codon}-{end_codon}"
        for ds in all_datasets:
            out_raw[ds] = str(sums[ds])

        collapsed.append(
            {
                "raw": out_raw,
                "gene": representative["gene"],
                "codon": representative["codon"],
            }
        )

    return collapsed


def collapse_rows(
    rows: List[Dict[str, object]],
    all_datasets: Sequence[str],
    selected_datasets: Sequence[str],
) -> List[Dict[str, object]]:
    # The collapsing logic works gene-by-gene so codons from different genes can
    # never be merged together accidentally.
    by_gene: Dict[str, List[Dict[str, object]]] = {}
    for row in rows:
        by_gene.setdefault(str(row["gene"]), []).append(row)

    out: List[Dict[str, object]] = []
    for gene in sorted(by_gene.keys()):
        out.extend(collapse_gene_rows(by_gene[gene], all_datasets, selected_datasets))

    return sorted(out, key=lambda r: (str(r["gene"]), int(r["codon"])))


def output_headers(headers: Sequence[str], selected_datasets: Sequence[str]) -> List[str]:
    # The output keeps all metadata columns through "Nucleotide Window", then
    # appends only the dataset columns that were actually selected for analysis.
    nuc_idx = headers.index(NUC_WINDOW_COL)
    metadata = list(headers[: nuc_idx + 1])
    if GENE_COL in metadata and PEAK_REGION_COL not in metadata:
        gene_idx = metadata.index(GENE_COL)
        metadata.insert(gene_idx + 1, PEAK_REGION_COL)
    return metadata + list(selected_datasets)


def main() -> None:
    # main() wires the whole workflow together:
    # - parse command-line options
    # - load the CSV and discover dataset columns
    # - filter qualifying rows
    # - collapse adjacent rows into regional peaks
    # - write the final CSV and print a short run summary
    args = parse_args()

    input_path = Path(args.input_csv)
    if not input_path.is_file():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    # If the user does not provide an output path, write beside the input file
    # with a "_collapsed" suffix so the original CSV is preserved.
    output_path = (
        Path(args.output_csv)
        if args.output_csv
        else input_path.with_name(f"{input_path.stem}_collapsed.csv")
    )

    include = split_csv_list(args.include_datasets)
    exclude = split_csv_list(args.exclude_datasets)

    with input_path.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        headers = list(reader.fieldnames or [])
        _nuc_idx, all_datasets = dataset_columns(headers)
        selected = select_datasets(all_datasets, include, exclude)

        # These validation checks catch contradictory command-line settings
        # early and explain the problem before any data is written.
        if not selected:
            raise ValueError("No datasets selected after include/exclude settings.")

        if args.min_datasets < 1:
            raise ValueError("--min-datasets must be >= 1")
        if args.min_datasets > len(selected):
            raise ValueError(
                f"--min-datasets ({args.min_datasets}) cannot exceed selected dataset count ({len(selected)})."
            )

        # First pass over the CSV: keep only codons that meet the user's
        # support threshold.
        kept = qualifying_rows(reader, all_datasets, selected, args.cutoff, args.min_datasets)

    # Second pass over the filtered rows: merge local codon neighborhoods into
    # broader regional peaks within each gene.
    collapsed = collapse_rows(kept, all_datasets, selected)
    out_headers = output_headers(headers, selected)

    with output_path.open("w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=out_headers, extrasaction="ignore")
        writer.writeheader()
        # Each collapsed item still contains the original row dictionary with
        # updated dataset sums, so writing is mostly a matter of streaming rows
        # back out in their final order.
        for row in collapsed:
            writer.writerow(row["raw"])

    # The printed summary helps the user verify how aggressive the filtering and
    # collapsing steps were for this run.
    print(f"Input rows qualifying cutoff: {len(kept)}")
    print(f"Output rows after collapsing: {len(collapsed)}")
    print(f"Selected datasets ({len(selected)}): {', '.join(selected)}")
    print(f"Wrote: {output_path}")


if __name__ == "__main__":
    main()
