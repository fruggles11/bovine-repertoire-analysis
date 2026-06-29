#!/usr/bin/env python3
"""
Extract ultra-long CDR3H3 sequences from bovine heavy chain majority consensus
FASTAs using direct J gene pattern matching.

Bypasses IgBLAST's CDR3 detection, which reports 'Total identifiable CDR3 = 0'
for ultra-long bovine CDR3H3 sequences because they exceed IgBLAST's heuristic
window. Instead:
  1. Uses IgBLAST's V gene call (which is reliable) from an existing AIRR TSV
  2. Directly searches for each IGHJ germline sequence in the consensus
  3. Identifies the conserved Trp codon (IMGT position 118) within the J gene
  4. Scans back for the conserved Cys codon (last codon of FR3) to define CDR3 start
  5. Translates and reports the full CDR3H3

Usage:
    custom_cdr3_filter.py \
        --airr_dir  ultralong_results/filtered \
        --fasta_dir results/5_majority_consensus \
        --germline  bovine_germline/Bos_taurus_IgHJ_gaps.fasta \
        --output_dir ultralong_custom_cdr3 \
        --min_cdr3_aa 40
"""

import argparse
import csv
import glob
import io
import re
import sys
from collections import Counter
from pathlib import Path


_ACGT = frozenset('ACGT')

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def revcomp(seq):
    comp = str.maketrans('ACGTacgtNnRrYyKkMmSsWwBbDdHhVv',
                         'TGCAtgcaNnYyRrMmKkSsWwVvHhDdBb')
    return seq.translate(comp)[::-1]


def translate(nt_seq):
    aa = []
    for i in range(0, len(nt_seq) - 2, 3):
        codon = nt_seq[i:i+3].upper()
        # Treat any codon with N as X
        if 'N' in codon or len(codon) < 3:
            aa.append('X')
        else:
            aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)


def parse_fasta(path):
    """Return dict of {seq_id: sequence}."""
    seqs = {}
    current_id = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        seqs[current_id] = ''.join(current_seq)
    return seqs


def load_j_germlines(fasta_path):
    """
    Load J gene sequences and record the position of the conserved Trp (TGG)
    within each gene.  Returns list of (name, clean_seq, trp_pos).
    """
    raw = parse_fasta(fasta_path)
    j_genes = []
    for header, seq in raw.items():
        name = header.split('|')[1] if '|' in header else header
        clean = re.sub(r'[.\-\s]', '', seq).upper()
        trp_pos = clean.find('TGG')
        if trp_pos == -1:
            continue  # skip J genes with no Trp (non-functional)
        j_genes.append((name, clean, trp_pos))
    return j_genes


def find_j_in_sequence(sequence, j_genes, anchor_len=20, min_identity=0.75):
    """
    Slide each J gene anchor across the sequence.
    Return (trp_position_in_sequence, j_name, identity) for the best match,
    or (None, None, 0) if nothing meets the threshold.
    """
    seq_u = sequence.upper()
    best = (None, None, 0.0)

    for j_name, j_seq, trp_pos in j_genes:
        # Anchor = first `anchor_len` nt of J gene starting from its Trp codon
        # (everything before the Trp is D-J junction noise)
        anchor_start = trp_pos
        anchor = j_seq[anchor_start: anchor_start + anchor_len]
        if len(anchor) < anchor_len:
            continue

        for i in range(len(seq_u) - len(anchor) + 1):
            window = seq_u[i: i + len(anchor)]
            # Only count unambiguous bases (A/C/G/T) toward identity;
            # all IUPAC ambiguity codes (N, M, Y, R, W, S, K, ...) are skipped
            matches = sum(
                1 for a, b in zip(anchor, window)
                if b in _ACGT and a == b
            )
            unambiguous = sum(1 for b in window if b in _ACGT)
            if unambiguous == 0:
                continue
            identity = matches / unambiguous
            if identity > best[2]:
                trp_in_seq = i  # anchor starts at the Trp codon
                best = (trp_in_seq, j_name, identity)

    if best[2] >= min_identity:
        return best
    return (None, None, 0.0)


def find_cys_before_trp(sequence, trp_pos, min_cdr3_aa=0, max_cdr3_nt=600):
    """
    Scan backwards from trp_pos in steps of 3 (same reading frame as Trp)
    looking for TGT or TGC (Cys codon).  Returns the position of the first
    Cys codon that gives CDR3H3 >= min_cdr3_aa, or None.

    Ultra-long bovine CDR3H3 regions contain internal Cys codons closer to the
    Trp (giving short apparent CDR3 lengths) that must be skipped.  Returning
    the first Cys that meets the length threshold gives the most conservative
    (shortest) CDR3H3 that still satisfies the filter.
    """
    seq_u = sequence.upper()
    for pos in range(trp_pos - 3, max(0, trp_pos - max_cdr3_nt) - 1, -3):
        codon = seq_u[pos: pos + 3]
        if codon in ('TGT', 'TGC'):
            cdr3_aa = (trp_pos + 3 - pos) // 3 - 2
            if cdr3_aa >= min_cdr3_aa:
                return pos
    return None


def extract_junction(sequence, cys_pos, trp_pos):
    """Extract nucleotide junction (Cys codon through Trp codon inclusive)."""
    return sequence[cys_pos: trp_pos + 3]


def cdrh3_len(junction_aa):
    """CDR3H3 length = junction_aa length minus Cys and Trp anchor residues."""
    return max(0, len(junction_aa) - 2)


def process_sequence(seq_id, sequence, j_genes,
                     min_cdr3_aa=40, min_j_identity=0.75):
    """
    Try both strands.  Return dict with CDR3 info or None if not found.
    """
    for strand, test_seq in [('+', sequence), ('-', revcomp(sequence))]:
        trp_pos, j_name, identity = find_j_in_sequence(
            test_seq, j_genes, min_identity=min_j_identity
        )
        if trp_pos is None:
            continue
        cys_pos = find_cys_before_trp(test_seq, trp_pos, min_cdr3_aa=min_cdr3_aa)
        if cys_pos is None:
            continue
        junction_nt = extract_junction(test_seq, cys_pos, trp_pos)
        junction_aa = translate(junction_nt)
        cdr3_len = cdrh3_len(junction_aa)
        if cdr3_len >= min_cdr3_aa:
            return {
                'strand':       strand,
                'j_call':       j_name,
                'j_identity':   round(identity, 4),
                'junction':     junction_nt,
                'junction_aa':  junction_aa,
                'junction_length': len(junction_nt),
                'cdrh3_length_aa': cdr3_len,
            }
    return None


def main():
    parser = argparse.ArgumentParser(
        description='Extract ultra-long CDR3H3 using direct J gene matching'
    )
    parser.add_argument('--airr_dir', required=True,
                        help='Directory with *_heavy_productive.tsv AIRR files '
                             '(from ultralong_cdrh3_filter.nf)')
    parser.add_argument('--fasta_dir', required=True,
                        help='Directory with *_heavy_majority.fasta files '
                             '(from bovine-igg-pipeline 5_majority_consensus/)')
    parser.add_argument('--germline', required=True,
                        help='Bovine IGHJ germline FASTA '
                             '(e.g. bovine_germline/Bos_taurus_IgHJ_gaps.fasta)')
    parser.add_argument('--output_dir', default='ultralong_custom_cdr3',
                        help='Output directory (default: ultralong_custom_cdr3)')
    parser.add_argument('--min_cdr3_aa', type=int, default=40,
                        help='Minimum CDR3H3 length in amino acids (default: 40)')
    parser.add_argument('--min_j_identity', type=float, default=0.75,
                        help='Minimum J gene search identity (default: 0.75)')
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading J gene germlines from {args.germline}")
    j_genes = load_j_germlines(args.germline)
    print(f"  Loaded {len(j_genes)} J genes with conserved Trp\n")

    # Load all majority consensus sequences keyed by barcode
    print(f"Loading majority consensus sequences from {args.fasta_dir}")
    fasta_seqs = {}  # barcode -> {seq_id -> sequence}
    for fasta_file in sorted(Path(args.fasta_dir).glob('*_heavy_majority.fasta')):
        barcode = fasta_file.stem.replace('_heavy_majority', '')
        fasta_seqs[barcode] = parse_fasta(fasta_file)
    print(f"  Loaded {len(fasta_seqs)} barcodes\n")

    # Load existing AIRR TSVs for V gene calls
    print(f"Loading AIRR annotations from {args.airr_dir}")
    airr_data = {}  # barcode -> row dict
    for tsv_file in sorted(Path(args.airr_dir).glob('*_heavy_productive.tsv')):
        barcode = tsv_file.stem.split('_heavy')[0]
        with open(tsv_file) as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        if rows:
            airr_data[barcode] = rows[0]
    print(f"  Loaded {len(airr_data)} AIRR annotations\n")

    # Process each barcode
    results = []
    for barcode in sorted(fasta_seqs.keys()):
        seqs = fasta_seqs[barcode]
        airr_row = airr_data.get(barcode, {})
        v_call = airr_row.get('v_call', '') or ''

        if v_call and 'IGHV' not in v_call:
            continue  # skip light chain sequences

        for seq_id, sequence in seqs.items():
            cdr3_info = process_sequence(
                seq_id, sequence, j_genes,
                min_cdr3_aa=args.min_cdr3_aa,
                min_j_identity=args.min_j_identity
            )
            if cdr3_info:
                results.append({
                    'barcode':        barcode,
                    'sequence_id':    seq_id,
                    'v_call':         v_call,
                    **cdr3_info,
                    'sequence':       sequence,
                })
                print(f"  {barcode}: CDR3H3 = {cdr3_info['cdrh3_length_aa']} aa  "
                      f"v={v_call}  j={cdr3_info['j_call']} "
                      f"(id={cdr3_info['j_identity']})  "
                      f"strand={cdr3_info['strand']}")

    print(f"\nTotal: {len(results)} sequences with CDR3H3 >= {args.min_cdr3_aa} aa")

    if not results:
        print("No sequences passed the filter.")
        sys.exit(0)

    results.sort(key=lambda r: r['cdrh3_length_aa'], reverse=True)

    # Write TSV
    tsv_out = out_dir / 'ultralong_cdrh3_custom.tsv'
    cols = ['barcode', 'sequence_id', 'cdrh3_length_aa', 'v_call', 'j_call',
            'j_identity', 'strand', 'junction_aa', 'junction_length',
            'junction', 'sequence']
    with open(tsv_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results)

    # Write FASTA
    fasta_out = out_dir / 'ultralong_cdrh3_custom.fasta'
    written = 0
    with open(fasta_out, 'w') as f:
        for r in results:
            seq = r['sequence']
            if r['strand'] == '-':
                seq = revcomp(seq)
            f.write(f">{r['barcode']}_{r['sequence_id']} "
                    f"v_call={r['v_call']} j_call={r['j_call']} "
                    f"cdrh3_aa={r['cdrh3_length_aa']} "
                    f"junction_aa={r['junction_aa']}\n")
            f.write(seq + '\n')
            written += 1

    print(f"\nOutput:")
    print(f"  {tsv_out} ({len(results)} entries)")
    print(f"  {fasta_out} ({written} sequences)")

    print(f"\nCDR3H3 length distribution:")
    dist = Counter(r['cdrh3_length_aa'] for r in results)
    for length in sorted(dist):
        bar = '#' * min(dist[length], 40)
        print(f"  {length:3d} aa: {bar} ({dist[length]})")


if __name__ == '__main__':
    main()
