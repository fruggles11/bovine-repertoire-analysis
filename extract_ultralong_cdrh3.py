#!/usr/bin/env python3
"""
Extract bovine sequences with ultra-long CDRH3 regions from AIRR repertoire TSV files.

Bovine ultra-long CDRH3 threshold: >= 50 aa (based on the bimodal distribution
observed in this dataset, with a dense cluster at 49-68 aa distinct from the
normal-length population peaking at 17-26 aa).
"""

import os
import csv
import glob
import argparse


def cdrh3_length_aa(junction_aa):
    """CDRH3 length = junction_aa length minus the two anchor residues (Cys + Trp/Phe)."""
    if not junction_aa:
        return None
    return len(junction_aa) - 2


def extract_ultralong(input_dir, output_file, threshold=50):
    pattern = os.path.join(input_dir, "barcode*_heavy_productive.tsv")
    tsv_files = sorted(glob.glob(pattern))

    if not tsv_files:
        print(f"No files matching pattern: {pattern}")
        return

    print(f"Found {len(tsv_files)} productive TSV files")
    print(f"CDRH3 length threshold: >= {threshold} aa\n")

    key_cols = [
        "sequence_id", "v_call", "d_call", "j_call", "c_call",
        "junction_aa", "junction_length", "junction",
        "np1_length", "np2_length",
        "stop_codon", "vj_in_frame", "productive", "locus",
        "sequence",
    ]

    total = 0
    rows_out = []
    barcode_counts = {}

    for tsv in tsv_files:
        barcode = os.path.basename(tsv).split("_heavy")[0]
        count = 0
        with open(tsv, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                jaa = row.get("junction_aa", "")
                length = cdrh3_length_aa(jaa)
                if length is not None and length >= threshold:
                    row["cdrh3_length_aa"] = length
                    row["barcode"] = barcode
                    rows_out.append(row)
                    count += 1
        barcode_counts[barcode] = count
        print(f"  {barcode}: {count} ultra-long CDRH3 sequences")
        total += count

    print(f"\nTotal: {total} sequences with CDRH3 >= {threshold} aa")

    if not rows_out:
        print("No sequences found above threshold.")
        return

    # Sort by CDRH3 length descending
    rows_out.sort(key=lambda r: r["cdrh3_length_aa"], reverse=True)

    # Write output
    out_cols = ["barcode", "sequence_id", "cdrh3_length_aa"] + key_cols
    # Only include cols that exist in the data
    available = set(rows_out[0].keys())
    out_cols = [c for c in out_cols if c in available]

    with open(output_file, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"\nOutput written to: {output_file}")

    # Print length distribution summary
    print("\nCDRH3 length distribution of extracted sequences:")
    from collections import Counter
    dist = Counter(r["cdrh3_length_aa"] for r in rows_out)
    for length in sorted(dist):
        bar = "#" * dist[length]
        print(f"  {length:3d} aa: {bar} ({dist[length]})")


def main():
    parser = argparse.ArgumentParser(description="Extract ultra-long CDRH3 bovine sequences")
    parser.add_argument(
        "--input-dir",
        default="repertoire_results/filtered",
        help="Directory containing barcode*_heavy_productive.tsv files"
    )
    parser.add_argument(
        "--output",
        default="ultralong_cdrh3_sequences.tsv",
        help="Output TSV file"
    )
    parser.add_argument(
        "--threshold",
        type=int,
        default=50,
        help="Minimum CDRH3 length in amino acids (default: 50)"
    )
    args = parser.parse_args()

    extract_ultralong(args.input_dir, args.output, args.threshold)


if __name__ == "__main__":
    main()
