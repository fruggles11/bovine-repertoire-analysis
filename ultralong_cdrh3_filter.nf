#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Ultra-Long CDR3H3 Extraction Pipeline
========================================================================================
    Identifies heavy chain sequences with ultra-long CDR3H3 regions from the
    majority consensus FASTAs produced by bovine-igg-pipeline.

    Steps:
      1. Build IgBLAST database from bovine germline FASTA files
      2. Annotate majority consensus FASTAs with V(D)J assignments via IgBLAST
      3. Convert to AIRR format with Change-O MakeDb
      4. Filter for CDR3H3 >= min_cdr3_aa amino acids
      5. Output a combined TSV and FASTA sorted by CDR3H3 length

    Usage:
      nextflow run ultralong_cdrh3_filter.nf \
          --fasta_dir /path/to/results/5_majority_consensus \
          --germline_db '/path/to/germlines/*.fasta'

    Options:
      --fasta_dir     Directory containing *_heavy_majority.fasta files
                      (default: ./results/5_majority_consensus)
      --germline_db   Glob pattern for bovine germline FASTA files from IMGT
                      (required — see README for download instructions)
      --min_cdr3_aa   Minimum CDR3H3 length in amino acids (default: 50)
      --results       Output directory (default: ./ultralong_results)
========================================================================================
*/

include { BUILD_IGBLAST_DB; IGBLAST_ANNOTATION; MAKEDB; FILTER_PRODUCTIVE } from './main.nf'

params.fasta_dir             = "${launchDir}/results/5_majority_consensus"
params.germline_db           = null
params.min_cdr3_aa           = 50
params.results               = "${launchDir}/ultralong_results"
params.skip_productive_filter = true   // productivity detection unreliable for bovine


// --------------------------------------------------------------- //
workflow {

    if ( !params.germline_db ) {
        error """
        ERROR: --germline_db is required.

        Download bovine germline sequences from IMGT/GENE-DB:
          1. Go to https://www.imgt.org/genedb/
          2. Select Species: "Bos taurus"
          3. Download IGHV, IGHD, IGHJ as FASTA files
          4. Save to a local directory (e.g., resources/germlines/)

        Then run with:
          --germline_db 'resources/germlines/*.fasta'
        """.stripIndent()
    }

    // Heavy chain majority consensus FASTAs from bovine-igg-pipeline
    ch_fasta = Channel
        .fromPath( "${params.fasta_dir}/*_heavy_majority.fasta" )
        .map { fasta ->
            def sample_id = fasta.baseName.replaceAll(/_majority$/, '')
            tuple( sample_id, fasta )
        }

    ch_germlines = Channel.fromPath( params.germline_db )

    // Build IgBLAST database
    BUILD_IGBLAST_DB( ch_germlines.collect() )

    // Annotate sequences with V(D)J gene assignments
    IGBLAST_ANNOTATION(
        ch_fasta,
        BUILD_IGBLAST_DB.out.db,
        BUILD_IGBLAST_DB.out.aux
    )

    // Convert IgBLAST output to AIRR-format TSV
    MAKEDB(
        IGBLAST_ANNOTATION.out,
        ch_germlines.collect()
    )

    // Pass-through productive filter (bovine productivity detection unreliable)
    FILTER_PRODUCTIVE( MAKEDB.out )

    // Collect all per-sample TSVs and extract ultra-long CDR3H3 sequences
    EXTRACT_ULTRALONG(
        FILTER_PRODUCTIVE.out
            .map { sample_id, tsv -> tsv }
            .collect()
    )

}
// --------------------------------------------------------------- //


// --------------------------------------------------------------- //
process EXTRACT_ULTRALONG {

    publishDir "${params.results}", mode: 'copy', overwrite: true

    input:
    path tsv_files

    output:
    path "ultralong_cdrh3.tsv"
    path "ultralong_cdrh3.fasta"

    script:
    """
    python3 << 'PYEOF'
import csv, glob, sys
from collections import Counter

threshold = ${params.min_cdr3_aa}
tsv_files = sorted(glob.glob('*_productive.tsv'))

print(f"Processing {len(tsv_files)} sample TSV files")
print(f"CDR3H3 threshold: >= {threshold} aa\\n")

rows_out = []
for tsv in tsv_files:
    barcode = tsv.split('_heavy')[0] if '_heavy' in tsv else tsv.replace('_productive.tsv', '')
    count = 0
    with open(tsv, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\\t')
        for row in reader:
            jaa = row.get('junction_aa', '') or ''
            # CDR3H3 length = junction_aa minus the two anchor residues (Cys + Trp/Phe)
            cdrh3_len = len(jaa) - 2 if len(jaa) >= 2 else 0
            if cdrh3_len >= threshold:
                row['cdrh3_length_aa'] = cdrh3_len
                row['barcode'] = barcode
                rows_out.append(row)
                count += 1
    print(f"  {barcode}: {count} sequences with CDR3H3 >= {threshold} aa")

print(f"\\nTotal: {len(rows_out)} ultra-long CDR3H3 sequences")

if not rows_out:
    print("No sequences passed the filter.")
    with open('ultralong_cdrh3.tsv', 'w') as f:
        f.write('sequence_id\\tv_call\\tjunction_aa\\tcdrh3_length_aa\\tsequence\\tbarcode\\n')
    open('ultralong_cdrh3.fasta', 'w').close()
    sys.exit(0)

# Sort by CDR3H3 length descending
rows_out.sort(key=lambda r: r['cdrh3_length_aa'], reverse=True)

# Write combined TSV (prioritised columns first, then any others)
priority_cols = ['barcode', 'sequence_id', 'cdrh3_length_aa', 'v_call', 'd_call', 'j_call',
                 'junction_aa', 'junction_length', 'junction', 'productive', 'sequence']
available = set(rows_out[0].keys())
out_cols = [c for c in priority_cols if c in available] + \
           [c for c in rows_out[0].keys() if c not in priority_cols]

with open('ultralong_cdrh3.tsv', 'w', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=out_cols, delimiter='\\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(rows_out)

# Write FASTA for import into Geneious
written = 0
with open('ultralong_cdrh3.fasta', 'w') as fh:
    for row in rows_out:
        seq = row.get('sequence', '') or ''
        if not seq:
            continue
        barcode  = row.get('barcode', '')
        seq_id   = row.get('sequence_id', '')
        v_call   = row.get('v_call', 'NA')
        cdr3_len = row.get('cdrh3_length_aa', 'NA')
        fh.write(f'>{barcode}_{seq_id} v_call={v_call} cdrh3_aa={cdr3_len}\\n')
        fh.write(seq + '\\n')
        written += 1

print(f"\\nOutput:")
print(f"  ultralong_cdrh3.tsv  ({len(rows_out)} entries)")
print(f"  ultralong_cdrh3.fasta ({written} sequences)")

# Print CDR3H3 length distribution
print("\\nCDR3H3 length distribution:")
dist = Counter(r['cdrh3_length_aa'] for r in rows_out)
for length in sorted(dist):
    bar = '#' * min(dist[length], 40)
    print(f"  {length:3d} aa: {bar} ({dist[length]})")
PYEOF
    """

}
// --------------------------------------------------------------- //
