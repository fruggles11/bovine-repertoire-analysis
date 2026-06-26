# Bovine IgG Repertoire Analysis Pipeline

Analyzes antibody repertoire diversity from bovine IgG sequences using the [Immcantation](https://immcantation.readthedocs.io/) framework. Includes a dedicated pipeline for extracting heavy chain sequences with ultra-long CDR3H3 regions, a hallmark feature of bovine immunoglobulins.

## Pipelines

| Pipeline | Script | Purpose |
|----------|--------|---------|
| Repertoire analysis | `main.nf` | Full V(D)J annotation, clonotype assignment, diversity analysis |
| Ultra-long CDR3H3 filter | `ultralong_cdrh3_filter.nf` | Extract heavy chains with CDR3H3 ≥ 50 aa for downstream use |

## Features

- **V(D)J annotation** using IgBLAST
- **Clonotype assignment** using Change-O
- **Diversity analysis** using Alakazam:
  - CDR3 length distribution
  - V and J gene usage
  - Clone size distribution
  - Hill diversity curves
  - Rarefaction curves
  - Simpson and Shannon diversity indices
- **Ultra-long CDR3H3 extraction** for identifying bovine-specific ultra-long CDR3H3 sequences

## Requirements

- Nextflow (≥21.04)
- Docker or Singularity

## Quick Start

```bash
# 1. Download bovine germlines from IMGT (see Germline Database section below)
#    Save FASTA files to a germlines/ directory

# 2. Run the pipeline on your bovine-igg-pipeline results
nextflow run fruggles11/bovine-repertoire-analysis \
    --germline_db '/path/to/germlines' \
    --input_dir '/path/to/results/4_consensus_sequences' \
    --results ./repertoire_results
```

The pipeline will automatically find all `*_unique.fasta` files recursively within the input directory.

## Input

The pipeline accepts FASTA files containing antibody sequences. These can be:
- `*_unique.fasta` files from the bovine-igg-pipeline
- `*_nogroup_unique.fasta` files
- Any FASTA with antibody V(D)J sequences

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_dir` | none | Directory to recursively search for `*_unique.fasta` files (recommended) |
| `--fasta_input` | `*_unique.fasta` | Glob pattern for input FASTA files (alternative to input_dir) |
| `--germline_db` | **required** | Path to bovine germline FASTAs (glob pattern or directory) |
| `--results` | `./results` | Output directory |
| `--clone_threshold` | 0.15 | Junction distance threshold for clonotype definition |

### Automatic Chain Filtering

The pipeline automatically detects the dominant chain type (heavy or light) for each barcode by counting sequences. The minority chain type is filtered out as cross-contamination.

For example, if barcode88 has 15,000 light chain sequences and 500 heavy chain sequences, the pipeline will:
- Keep: `barcode88_light` (dominant)
- Filter out: `barcode88_heavy` (contamination)

This works automatically with any number of barcodes - no manual configuration needed.

## Output

```
results/
├── igblast/           # Raw IgBLAST annotations
├── airr/              # AIRR-formatted sequence data
├── filtered/          # Productive sequences only
├── clones/            # Clone assignments
├── diversity/
│   ├── plots/         # PDF and PNG visualizations
│   │   ├── cdr3_length_distribution.png
│   │   ├── v_gene_usage.png
│   │   ├── j_gene_usage.png
│   │   ├── clone_size_distribution.png
│   │   ├── diversity_curve.png
│   │   └── rarefaction_curve.png
│   └── stats/         # CSV statistics
│       ├── basic_stats.csv
│       ├── v_gene_usage.csv
│       ├── j_gene_usage.csv
│       ├── clone_sizes.csv
│       └── diversity_summary.csv
└── repertoire_report.html
```

## Germline Database (Required)

The pipeline requires bovine immunoglobulin germline sequences. **Manual download is required** because Immcantation's `fetch_imgtdb.sh` does not support bovine.

> **⚠️ IMPORTANT: You must download IMGT-gapped sequences.** Ungapped sequences will not work correctly - CDR3/junction regions will not be identified and all downstream analysis will fail.

### Download from IMGT/GENE-DB

1. Go to [IMGT/GENE-DB](https://www.imgt.org/genedb/)
2. Select **Species**: "Bos taurus"
3. Select **Molecular component**: "IG" (immunoglobulin)
4. Select **Functionality**: "F" (functional) or leave as default
5. **CRITICAL**: Check the box for **"With IMGT gaps"** before downloading
6. For each gene type, click "FASTA" to download:

   **Heavy chain (required):**
   - **IGHV** - Variable genes
   - **IGHD** - Diversity genes
   - **IGHJ** - Joining genes

   **Light chain (if analyzing light chains):**
   - **IGKV**, **IGKJ** - Kappa chain
   - **IGLV**, **IGLJ** - Lambda chain

7. Save all FASTA files to a single directory (e.g., `germlines/`)
8. Provide the path via `--germline_db './germlines/*.fasta'` or simply `--germline_db './germlines'`

### Why IMGT-gapped sequences?

IMGT gaps are standardized insertions (shown as dots `.` in the sequence) that maintain consistent numbering across all immunoglobulin sequences. These gaps are required for:
- Proper CDR3/junction identification
- Accurate V(D)J boundary detection
- Correct reading frame determination

Without gaps, the pipeline will run but produce empty CDR3 data and incorrect diversity metrics.

### Example

```bash
mkdir -p germlines
# Download IMGT-gapped FASTA files and save to germlines/
# Files should look like: Bos_taurus_IGHV_gapped.fasta, etc.
# Then run:
nextflow run fruggles11/bovine-repertoire-analysis \
    --germline_db './germlines/*.fasta' \
    --input_dir '/path/to/results/4_consensus_sequences'
```

---

## Ultra-Long CDR3H3 Filter (`ultralong_cdrh3_filter.nf`)

Bovine heavy chains uniquely produce CDR3H3 regions up to 70+ amino acids long — far exceeding the ~15 aa typical in humans. This pipeline identifies these ultra-long sequences from the majority consensus FASTAs output by [bovine-igg-pipeline](https://github.com/fruggles11/bovine-igg-pipeline).

Only sequences with a confirmed IGHV gene call are reported — IgBLAST occasionally assigns light chain germlines (IGLV) to heavy chain sequences, so this filter prevents false positives.

### Quick Start

```bash
nextflow run fruggles11/bovine-repertoire-analysis/ultralong_cdrh3_filter.nf \
    --fasta_dir /path/to/results/5_majority_consensus \
    --germline_db '/path/to/germlines/*.fasta'
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fasta_dir` | `./results/5_majority_consensus` | Directory containing `*_heavy_majority.fasta` files from bovine-igg-pipeline |
| `--germline_db` | **required** | Glob pattern for bovine germline FASTAs from IMGT (same files as main pipeline) |
| `--min_cdr3_aa` | `40` | Minimum CDR3H3 length in amino acids |
| `--results` | `./ultralong_results` | Output directory |

The 40 aa default threshold is based on the bimodal distribution observed in bovine heavy chain repertoire data, where the ultra-long population is clearly separated from the normal-length population (which peaks at 17–26 aa).

### Output

```
ultralong_results/
├── ultralong_cdrh3.fasta   # Filtered sequences sorted by CDR3H3 length (for Geneious import)
├── ultralong_cdrh3.tsv     # Full AIRR annotation table with CDR3H3 length column
├── igblast/                # Per-sample IgBLAST annotations
├── airr/                   # Per-sample AIRR-format TSVs
└── filtered/               # Per-sample productive TSVs
```

The FASTA headers include V gene assignment and CDR3H3 length:
```
>barcode41_barcode41 v_call=IGHV1-7*02 cdrh3_aa=62
```

### Notes on CDR3H3 detection

IgBLAST reports `Total identifiable CDR3 = 0` for some ultra-long sequences because the CDR3 length exceeds its internal heuristic window. In practice this means sequences with very long CDR3s (>55 aa) are reliably detected, while sequences in the 40–50 aa range may occasionally be missed or have slightly underestimated lengths. Running with `--min_cdr3_aa 40` captures the full ultra-long population with minimal false positives when the IGHV filter is applied.

### Utilities

**`extract_ultralong_cdrh3.py`** — standalone script for filtering an existing directory of `*_heavy_productive.tsv` files without re-running IgBLAST:

```bash
python3 extract_ultralong_cdrh3.py \
    --input-dir ultralong_results/filtered \
    --output ultralong_cdrh3_sequences.tsv \
    --threshold 40
```

---

## Diversity Metrics

The pipeline calculates:

- **Species richness**: Number of unique clonotypes
- **Shannon index**: Diversity accounting for abundance
- **Simpson index**: Probability two random sequences are different clones
- **Chao1**: Estimated true richness including unseen clones
- **Hill numbers**: Unified diversity framework at different q values

## Citation

If you use this pipeline, please cite:

- **Immcantation**: Gupta NT, et al. (2015) Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data. Bioinformatics.
- **Alakazam**: Stern JNH, et al. (2014) B cells populating the multiple sclerosis brain mature in the draining cervical lymph nodes. Science Translational Medicine.

## Tools Used

- [IgBLAST](https://www.ncbi.nlm.nih.gov/igblast/) - V(D)J gene assignment
- [Change-O](https://changeo.readthedocs.io/) - Clonal assignment
- [Alakazam](https://alakazam.readthedocs.io/) - Diversity analysis
- [Shazam](https://shazam.readthedocs.io/) - Statistical analysis
