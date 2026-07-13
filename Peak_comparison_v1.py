"""
Peak comparison utilities.

This module loads multiple ribosome-density pickle files, normalizes codon
densities per gene per file, and reports codon positions where at least one
dataset exceeds a user-provided enrichment cutoff.

Primary entrypoint:
    peak_comparison_wrapper(txtinfile)

Example command-line usage:
    python3 /path/to/mammaliancode/peakcomparisoncaller.py

This uses the settings stored in:
    /path/to/mammaliancode/inputfiles/peakcomparison.txt

How this script runs, at a high level:
1) Read an INPUT_* text file to collect file paths and runtime settings.
2) Load each pkl.gzip footprint dictionary into memory.
3) For each gene, shift the 5' read-end track by the requested offset.
4) Convert the shifted nucleotide track into codon-level ORF densities.
5) Normalize each gene within each file by that gene's mean non-zero codon density.
6) Write one CSV row for any codon where at least one dataset exceeds cutoff.
"""

import csv
import gzip
import os
import pickle
import sys
import time
from concurrent.futures import ThreadPoolExecutor

try:
    import resource
except Exception:
    resource = None


# Standard DNA codon table for short translated display windows in output.
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_nt(nt_seq):
    """Translate a nucleotide string in-frame (3 nt blocks) for display."""
    # The output CSV includes a short amino-acid window to help humans identify
    # the local context around a reported codon peak.
    nt_seq = nt_seq.upper()
    aa = []
    for i in range(0, len(nt_seq) - 2, 3):
        codon = nt_seq[i:i + 3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)


def get_rss_mb():
    """
    Return peak resident memory usage in MB, if available on this platform.
    """
    if resource is None:
        return None
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return float(rss) / (1024.0 * 1024.0)
    return float(rss) / 1024.0


def shift_counts_5prime(counts, offset):
    """Shift 5' end-aligned counts by offset with zero padding at boundaries."""
    # The pkl files store 5' end-aligned counts. Peak calling is done on a
    # site-assigned profile, so this shift moves each read to the requested
    # inferred ribosome position before ORF slicing begins.
    shifted = [0.0 for _ in range(len(counts))]
    for idx, value in enumerate(counts):
        new_idx = idx + offset
        if 0 <= new_idx < len(counts):
            shifted[new_idx] += float(value)
    return shifted


def load_single_peak_file(file_name, timer_enabled=False):
    """
    Load one .pkl.gzip footprint file and return:
        (file_name, footprints_dict, elapsed_seconds, compressed_size_mb)
    """
    if not os.path.isfile(file_name):
        print("File '{}' not found.".format(file_name))
        return None

    start_time = time.time()
    try:
        # Each input file is a gzip-compressed pickle containing one dictionary
        # keyed by transcript/gene identifier.
        with gzip.open(file_name, "rb") as loaded_data:
            footprints = pickle.load(loaded_data)

        elapsed = time.time() - start_time
        file_size_mb = float(os.path.getsize(file_name)) / (1024.0 * 1024.0)

        if timer_enabled:
            print(
                "Loaded {} in {:.2f}s ({:.1f} MB compressed).".format(
                    file_name, elapsed, file_size_mb
                )
            )
        else:
            print("Successfully loaded: {}".format(file_name))
        return (file_name, footprints, elapsed, file_size_mb)
    except Exception as e:
        print("Error processing '{}': {}".format(file_name, e))
        return None


def load_peak_files(file_list, load_workers=1, timer_enabled=False, memory_enabled=False):
    """
    Load all input files.

    If load_workers > 1, loading is done in parallel threads.
    Returns:
        results[file_name] = footprints_dict
        loaded_files = list of successfully loaded files
    """
    results = {}
    loaded_files = []
    load_start = time.time()

    # Loading can dominate runtime for large footprint dictionaries, so the
    # script optionally parallelizes file IO/decompression across threads.
    if load_workers <= 1:
        for file_name in file_list:
            item = load_single_peak_file(file_name, timer_enabled=timer_enabled)
            if item is None:
                continue
            results[item[0]] = item[1]
            loaded_files.append(item[0])
    else:
        with ThreadPoolExecutor(max_workers=load_workers) as executor:
            for item in executor.map(load_single_peak_file, file_list, [timer_enabled] * len(file_list)):
                if item is None:
                    continue
                results[item[0]] = item[1]
                loaded_files.append(item[0])

    if timer_enabled:
        print(
            "Total load time: {:.2f}s for {} files.".format(
                time.time() - load_start, len(loaded_files)
            )
        )
    if memory_enabled:
        rss_mb = get_rss_mb()
        if rss_mb is not None:
            print("Peak RSS after loading: {:.1f} MB".format(rss_mb))

    return results, loaded_files


def resolve_output_filename(outfile):
    """
    Resolve output csv path strictly from INPUT_outfile.
    """
    if outfile is None or outfile.strip() == "":
        raise ValueError("INPUT_outfile is required in peakcomparison.txt and cannot be blank.")

    output_name = outfile.strip()
    if not output_name.endswith(".csv"):
        output_name += ".csv"
    print("Output file is: " + output_name)
    return output_name


def find_significant_values(
    results,
    files,
    offset,
    cutoff,
    outfile=None,
    timer_enabled=False,
    memory_enabled=False,
    progress_interval=2000,
):
    """
    Core comparison routine.

    Steps per gene:
    1) Validate gene exists and CDS lengths match across all datasets.
    2) Build per-file codon densities and per-file non-zero codon mean.
    3) Emit rows for codons where any dataset exceeds `cutoff`.
    """
    active_files = [f for f in files if f in results]
    if len(active_files) == 0:
        print("No valid input files were loaded. Exiting.")
        return

    analysis_start = time.time()
    genes_seen = 0
    output_rows = 0

    output_filename = resolve_output_filename(outfile)
    with open(output_filename, "w", newline="") as file_handle:
        writer = csv.writer(file_handle)

        # Header: metadata columns + one normalized-density column per dataset.
        header = ["Gene", "Codon Number", "Gene_id", "Codon window", "Nucleotide Window"]
        for file_name in active_files:
            header.append(os.path.basename(file_name).replace(".pkl.gzip", ""))
        writer.writerow(header)

        base_results = results[active_files[0]]
        file_count = len(active_files)

        # The first loaded file acts as the reference list of genes and CDS
        # bounds. Every other file must contain the same gene with the same ORF
        # length or that gene is skipped for this comparison.
        for gene_key in base_results:
            genes_seen += 1
            if timer_enabled and progress_interval > 0 and genes_seen % progress_interval == 0:
                print("Processed {} genes in {:.2f}s".format(genes_seen, time.time() - analysis_start))

            length_match = True
            gene_reads = False

            # Use base file gene bounds for sequence context and target codon count.
            base_start = int(base_results[gene_key][3])
            base_stop = int(base_results[gene_key][4])
            nt_length = base_stop - base_start
            aa_length = int(nt_length / 3)

            nt_orf_densities_by_file = [None] * file_count
            average_densities_by_file = [1.0] * file_count

            # First pass for this gene: build the site-shifted ORF density track
            # for each file and compute the per-gene normalization denominator.
            #
            # The denominator is the mean of non-zero codon sums within that
            # gene/file. This makes the later codon values relative enrichment
            # scores rather than raw counts.
            for idx, file_name in enumerate(active_files):
                if gene_key not in results[file_name]:
                    length_match = False
                    break

                gene_data = results[file_name][gene_key]
                file_start = int(gene_data[3])
                file_stop = int(gene_data[4])
                file_nt_length = file_stop - file_start
                if file_nt_length != nt_length:
                    length_match = False
                    break

                # Shift the full nucleotide track before slicing the ORF. This
                # avoids the negative-slice bug that can drop genes whose CDS
                # starts near the 5' end of the transcript.
                shifted_all_5 = shift_counts_5prime(gene_data[2]["all_5"], offset)
                nt_orf_density = shifted_all_5[file_start:file_stop]
                if len(nt_orf_density) != file_nt_length:
                    length_match = False
                    break

                nt_orf_densities_by_file[idx] = nt_orf_density

                # Mean of non-zero codon sums (normalization denominator).
                # Genes with no ORF reads in every file are skipped later.
                non_zero_sum = 0.0
                non_zero_count = 0
                for aa_idx in range(aa_length):
                    density = sum(nt_orf_density[aa_idx * 3:aa_idx * 3 + 3])
                    if density != 0.0:
                        non_zero_sum += density
                        non_zero_count += 1

                if non_zero_count == 0:
                    avg_density = 1.0
                else:
                    avg_density = non_zero_sum / non_zero_count
                    gene_reads = True

                average_densities_by_file[idx] = avg_density

            if (not length_match) or (not gene_reads):
                continue

            # Second pass for this gene: test each codon and emit a row if at
            # least one dataset crosses the enrichment cutoff.
            base_nt_seq = base_results[gene_key][1]
            base_alias = base_results[gene_key][0]
            for aa_idx in range(aa_length):
                codon_above_cutoff = False
                codon_values = []

                for idx in range(file_count):
                    # Codon density is the sum of the three nucleotides that
                    # belong to the codon after site assignment. Dividing by the
                    # per-gene mean produces the reported enrichment value.
                    codon_density = sum(
                        nt_orf_densities_by_file[idx][aa_idx * 3:aa_idx * 3 + 3]
                    )
                    codon_value = codon_density / average_densities_by_file[idx]
                    codon_values.append(codon_value)
                    if codon_value > cutoff:
                        codon_above_cutoff = True

                if not codon_above_cutoff:
                    continue

                # The output keeps a short local sequence window so downstream
                # review can connect a numeric peak with its sequence context.
                # This 21-nt window is legacy behavior retained for compatibility.
                nt_seq = base_nt_seq[aa_idx * 3 + base_start - 15:aa_idx * 3 + 6 + base_start]
                output = [base_alias, aa_idx + 1, gene_key, translate_nt(nt_seq), nt_seq]
                output.extend(codon_values)
                writer.writerow(output)
                output_rows += 1

    if timer_enabled:
        print(
            "Analysis time: {:.2f}s | genes processed: {} | output rows: {}".format(
                time.time() - analysis_start, genes_seen, output_rows
            )
        )
    if memory_enabled:
        rss_mb = get_rss_mb()
        if rss_mb is not None:
            print("Peak RSS after analysis: {:.1f} MB".format(rss_mb))


def parse_input_file(txtinfile):
    """
    Parse INPUT_* style config file used across this repository.
    """
    # The repository uses a simple alternating format:
    # INPUT_name
    # value
    #
    # Blank lines and comment lines are ignored.
    with open(txtinfile) as handle:
        input_params = [line.strip() for line in handle.readlines()]

    variables = {}
    nextlineisinput = 0
    for inputline in input_params:
        if inputline == "" or inputline[0] == "#":
            continue

        if nextlineisinput != 0:
            variables[nextlineisinput] = inputline
            nextlineisinput = 0

        if inputline[0:5] == "INPUT":
            nextlineisinput = inputline.split("_")[-1]

    return variables


def peak_comparison_wrapper(txtinfile):
    """
    Wrapper entrypoint:
    1) Parse txt input file,
    2) Load all data files,
    3) Run comparison and write csv.
    """
    total_start = time.time()
    variables = parse_input_file(txtinfile)

    # Echo the parsed inputs up front so the run log captures exactly which
    # files and thresholds produced the output CSV.
    print("Variables collected:")
    for key in variables.keys():
        print(key)
        print(variables[key])
        print("")

    file_list = [x.strip() for x in variables["picklenames"].split(",") if x.strip() != ""]
    offset = int(variables.get("offset", "43"))
    cutoff = float(variables.get("cutoff", "5"))
    outfile = variables.get("outfile", "")
    load_workers = int(variables.get("loadworkers", "1"))
    timer_enabled = int(variables.get("timer", "1")) == 1
    memory_enabled = int(variables.get("memoryreport", "1")) == 1
    progress_interval = int(variables.get("progressinterval", "2000"))

    if timer_enabled:
        print("Load workers: {}".format(load_workers))

    # The wrapper is intentionally small: load data once, run the comparison
    # once, and let the lower-level functions handle all filtering/output logic.
    results, loaded_files = load_peak_files(
        file_list,
        load_workers=load_workers,
        timer_enabled=timer_enabled,
        memory_enabled=memory_enabled,
    )

    find_significant_values(
        results,
        loaded_files,
        offset,
        cutoff,
        outfile=outfile,
        timer_enabled=timer_enabled,
        memory_enabled=memory_enabled,
        progress_interval=progress_interval,
    )

    if timer_enabled:
        print("Total runtime: {:.2f}s".format(time.time() - total_start))
