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

### Download from IMGT/GENE-DB

1. Go to [IMGT/GENE-DB](https://www.imgt.org/genedb/)
2. Select **Species**: "Bos taurus"
3. Select **Molecular component**: "IG" (immunoglobulin)
4. For each gene type, download the nucleotide sequences:
   - **IGHV** (heavy chain variable)
   - **IGHD** (heavy chain diversity)
   - **IGHJ** (heavy chain joining)
   - Optionally: IGKV, IGKJ, IGLV, IGLJ for light chains
5. Save all FASTA files to a single directory (e.g., `germlines/`)
6. Provide the path via `--germline_db './germlines/*.fasta'`

### Example

```bash
mkdir -p germlines
# Download FASTA files from IMGT and save to germlines/
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
