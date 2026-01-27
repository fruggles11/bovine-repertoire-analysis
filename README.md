# Bovine IgG Repertoire Analysis Pipeline

Analyzes antibody repertoire diversity from bovine IgG sequences using the [Immcantation](https://immcantation.readthedocs.io/) framework.

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

## Requirements

- Nextflow (≥21.04)
- Docker or Singularity

## Quick Start

```bash
# 1. Download bovine germlines from IMGT (see Germline Database section below)
#    Save FASTA files to a germlines/ directory

# 2. From the bovine-igg-pipeline results directory
cd /path/to/results/4_consensus_sequences

# 3. Collect all unique FASTA files
rm -rf analysis_input && mkdir -p analysis_input && find . -path ./analysis_input -prune -o -name "*_unique.fasta" -exec cp {} ./analysis_input/ \;

# 4. Run the pipeline (germline_db is required for bovine)
nextflow run fruggles11/bovine-repertoire-analysis \
    --germline_db '/path/to/germlines/*.fasta' \
    --fasta_input 'analysis_input/*_unique.fasta' \
    --results ./repertoire_results
```

## Input

The pipeline accepts FASTA files containing antibody sequences. These can be:
- `*_unique.fasta` files from the bovine-igg-pipeline
- `*_nogroup_unique.fasta` files
- Any FASTA with antibody V(D)J sequences

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fasta_input` | `*_unique.fasta` | Glob pattern for input FASTA files |
| `--germline_db` | **required** | Path to bovine germline FASTAs (glob pattern) |
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
8. Provide the path via `--germline_db './germlines/*.fasta'`

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
    --fasta_input 'analysis_input/*_unique.fasta'
```

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
