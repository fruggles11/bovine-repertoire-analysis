#!/usr/bin/env python3
"""
Find ultra-long CDR3H3 sequences from bovine heavy chain data.

Combines two detection methods and merges the results:

  Method 1 — IgBLAST: reads CDR3H3 length directly from IgBLAST's
    junction_aa annotation in AIRR TSV files. Reliable when IgBLAST can
    annotate the CDR3, but silently returns empty junction_aa for sequences
    where the CDR3 exceeds its internal heuristic window (~55 aa).

  Method 2 — J gene anchor: searches the majority consensus FASTA directly
    for the conserved Trp codon (IMGT position 118) via J germline matching,
    then scans back for the nearest Cys codon. Catches additional sequences
    IgBLAST misses, subject to the limitation that ultra-long CDR3H3 regions
    often contain internal Cys codons — the closest Cys to the Trp may be an
    internal residue rather than the FR3-boundary Cys, causing some genuine
    ultra-long CDR3s to be reported with underestimated lengths.

  Merge strategy: IgBLAST results take precedence (more accurate lengths).
    Method 2 is run only on barcodes IgBLAST did not detect.

Usage:
    find_ultralong_cdrh3.py \
        --airr_dir  ultralong_results/filtered \
        --fasta_dir results/5_majority_consensus \
        --germline  bovine_germline/Bos_taurus_IgHJ_gaps.fasta \
        --output_dir ultralong_cdrh3 \
        --min_cdr3_aa 40
"""

import argparse
import csv
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
        if 'N' in codon or len(codon) < 3:
            aa.append('X')
        else:
            aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)


def parse_fasta(path):
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


# ── Method 1: IgBLAST-based detection ────────────────────────────────────────

def detect_igblast(airr_dir, min_cdr3_aa):
    """
    Read *_heavy_productive.tsv files and return rows where IgBLAST annotated
    a CDR3H3 >= min_cdr3_aa with an IGHV v_call.
    Returns dict: barcode -> unified result dict.
    """
    results = {}
    for tsv_file in sorted(Path(airr_dir).glob('*_heavy_productive.tsv')):
        barcode = tsv_file.stem.split('_heavy')[0]
        with open(tsv_file) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                if 'IGHV' not in (row.get('v_call', '') or ''):
                    continue
                jaa = row.get('junction_aa', '') or ''
                length = max(0, len(jaa) - 2) if jaa else None
                if length is not None and length >= min_cdr3_aa:
                    results[barcode] = {
                        'barcode':         barcode,
                        'sequence_id':     row.get('sequence_id', barcode),
                        'cdrh3_length_aa': length,
                        'source':          'igblast',
                        'v_call':          row.get('v_call', ''),
                        'd_call':          row.get('d_call', ''),
                        'j_call':          row.get('j_call', ''),
                        'junction_aa':     jaa,
                        'junction_length': int(row.get('junction_length', 0) or 0),
                        'junction':        row.get('junction', ''),
                        'sequence':        row.get('sequence', ''),
                    }
    return results


# ── Method 2: J gene anchor detection ────────────────────────────────────────

def load_j_germlines(fasta_path):
    raw = parse_fasta(fasta_path)
    j_genes = []
    for header, seq in raw.items():
        name = header.split('|')[1] if '|' in header else header
        clean = re.sub(r'[.\-\s]', '', seq).upper()
        trp_pos = clean.find('TGG')
        if trp_pos == -1:
            continue
        j_genes.append((name, clean, trp_pos))
    return j_genes


def find_j_in_sequence(sequence, j_genes, anchor_len=20, min_identity=0.75):
    seq_u = sequence.upper()
    best = (None, None, 0.0)
    for j_name, j_seq, trp_pos in j_genes:
        anchor = j_seq[trp_pos: trp_pos + anchor_len]
        if len(anchor) < anchor_len:
            continue
        for i in range(len(seq_u) - len(anchor) + 1):
            window = seq_u[i: i + len(anchor)]
            matches = sum(1 for a, b in zip(anchor, window) if b in _ACGT and a == b)
            unambiguous = sum(1 for b in window if b in _ACGT)
            if unambiguous == 0:
                continue
            identity = matches / unambiguous
            if identity > best[2]:
                best = (i, j_name, identity)
    if best[2] >= min_identity:
        return best
    return (None, None, 0.0)


def find_cys_before_trp(sequence, trp_pos, max_cdr3_nt=600):
    seq_u = sequence.upper()
    for pos in range(trp_pos - 3, max(0, trp_pos - max_cdr3_nt) - 1, -3):
        codon = seq_u[pos: pos + 3]
        if codon in ('TGT', 'TGC'):
            return pos
    return None


def detect_anchor(fasta_dir, airr_dir, j_genes, skip_barcodes,
                  min_cdr3_aa, min_j_identity):
    """
    Search consensus FASTAs for CDR3H3 via J gene anchor matching.
    Only processes barcodes not already in skip_barcodes (IgBLAST hits).
    Returns dict: barcode -> unified result dict.
    """
    # Load consensus FASTAs
    fasta_seqs = {}
    for fasta_file in sorted(Path(fasta_dir).glob('*_heavy_majority.fasta')):
        barcode = fasta_file.stem.replace('_heavy_majority', '')
        if barcode not in skip_barcodes:
            fasta_seqs[barcode] = parse_fasta(fasta_file)

    # Load v_calls from AIRR TSVs (first row per barcode)
    v_calls = {}
    for tsv_file in sorted(Path(airr_dir).glob('*_heavy_productive.tsv')):
        barcode = tsv_file.stem.split('_heavy')[0]
        if barcode in skip_barcodes:
            continue
        with open(tsv_file) as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        if rows:
            v_calls[barcode] = rows[0].get('v_call', '') or ''

    results = {}
    for barcode, seqs in sorted(fasta_seqs.items()):
        v_call = v_calls.get(barcode, '')
        if v_call and 'IGHV' not in v_call:
            continue

        for seq_id, sequence in seqs.items():
            for strand, test_seq in [('+', sequence), ('-', revcomp(sequence))]:
                trp_pos, j_name, identity = find_j_in_sequence(
                    test_seq, j_genes, min_identity=min_j_identity
                )
                if trp_pos is None:
                    continue
                cys_pos = find_cys_before_trp(test_seq, trp_pos)
                if cys_pos is None:
                    continue
                junction_nt = test_seq[cys_pos: trp_pos + 3]
                junction_aa = translate(junction_nt)
                cdr3_len = max(0, len(junction_aa) - 2)
                if cdr3_len >= min_cdr3_aa:
                    oriented_seq = sequence if strand == '+' else revcomp(sequence)
                    results[barcode] = {
                        'barcode':         barcode,
                        'sequence_id':     seq_id,
                        'cdrh3_length_aa': cdr3_len,
                        'source':          'anchor',
                        'v_call':          v_call,
                        'd_call':          '',
                        'j_call':          j_name,
                        'junction_aa':     junction_aa,
                        'junction_length': len(junction_nt),
                        'junction':        junction_nt,
                        'sequence':        oriented_seq,
                    }
                    break
            if barcode in results:
                break

    return results


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Find ultra-long CDR3H3 sequences using IgBLAST annotations '
                    'and J gene anchor matching'
    )
    parser.add_argument('--airr_dir', required=True,
                        help='Directory with *_heavy_productive.tsv AIRR files '
                             '(from ultralong_cdrh3_filter.nf filtered/ output)')
    parser.add_argument('--fasta_dir', required=True,
                        help='Directory with *_heavy_majority.fasta files '
                             '(from bovine-igg-pipeline 5_majority_consensus/)')
    parser.add_argument('--germline', required=True,
                        help='Bovine IGHJ germline FASTA '
                             '(e.g. bovine_germline/Bos_taurus_IgHJ_gaps.fasta)')
    parser.add_argument('--output_dir', default='ultralong_cdrh3',
                        help='Output directory (default: ultralong_cdrh3)')
    parser.add_argument('--min_cdr3_aa', type=int, default=40,
                        help='Minimum CDR3H3 length in amino acids (default: 40)')
    parser.add_argument('--min_j_identity', type=float, default=0.75,
                        help='Minimum J gene anchor match identity (default: 0.75)')
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Method 1: IgBLAST
    print('Method 1: IgBLAST annotation')
    igblast_hits = detect_igblast(args.airr_dir, args.min_cdr3_aa)
    for bc, r in sorted(igblast_hits.items(), key=lambda x: -x[1]['cdrh3_length_aa']):
        print(f'  {bc}: {r["cdrh3_length_aa"]} aa  v={r["v_call"]}  j={r["j_call"]}')
    print(f'  → {len(igblast_hits)} sequences\n')

    # Method 2: J gene anchor (only on barcodes IgBLAST missed)
    print('Method 2: J gene anchor matching')
    print(f'  Loading J germlines from {args.germline}')
    j_genes = load_j_germlines(args.germline)
    print(f'  Loaded {len(j_genes)} J genes')
    anchor_hits = detect_anchor(
        args.fasta_dir, args.airr_dir, j_genes,
        skip_barcodes=set(igblast_hits),
        min_cdr3_aa=args.min_cdr3_aa,
        min_j_identity=args.min_j_identity,
    )
    for bc, r in sorted(anchor_hits.items(), key=lambda x: -x[1]['cdrh3_length_aa']):
        print(f'  {bc}: {r["cdrh3_length_aa"]} aa  v={r["v_call"]}  j={r["j_call"]}')
    print(f'  → {len(anchor_hits)} sequences\n')

    # Merge: IgBLAST preferred, anchor as fallback
    merged = {**anchor_hits, **igblast_hits}
    results = sorted(merged.values(), key=lambda r: r['cdrh3_length_aa'], reverse=True)

    print(f'Total: {len(results)} sequences with CDR3H3 >= {args.min_cdr3_aa} aa')
    print(f'  {len(igblast_hits)} from IgBLAST  |  {len(anchor_hits)} from J anchor\n')

    if not results:
        print('No sequences passed the filter.')
        sys.exit(0)

    # Write TSV
    cols = ['barcode', 'sequence_id', 'cdrh3_length_aa', 'source',
            'v_call', 'd_call', 'j_call', 'junction_aa', 'junction_length',
            'junction', 'sequence']
    tsv_out = out_dir / 'ultralong_cdrh3.tsv'
    with open(tsv_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results)

    # Write FASTA
    fasta_out = out_dir / 'ultralong_cdrh3.fasta'
    with open(fasta_out, 'w') as f:
        for r in results:
            f.write(f">{r['barcode']} v_call={r['v_call']} "
                    f"cdrh3_aa={r['cdrh3_length_aa']} source={r['source']}\n")
            f.write(r['sequence'] + '\n')

    print(f'Output:')
    print(f'  {tsv_out}')
    print(f'  {fasta_out}')

    print(f'\nCDR3H3 length distribution:')
    dist = Counter(r['cdrh3_length_aa'] for r in results)
    for length in sorted(dist):
        bar = '#' * min(dist[length], 40)
        print(f'  {length:3d} aa: {bar} ({dist[length]})')


if __name__ == '__main__':
    main()
