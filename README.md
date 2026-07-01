# Bovine IgG Repertoire Analysis Pipeline

Analyzes antibody repertoire diversity from bovine IgG sequences using the [Immcantation](https://immcantation.readthedocs.io/) framework. Includes a dedicated pipeline for extracting heavy chain sequences with ultra-long CDR3H3 regions, a hallmark feature of bovine immunoglobulins.

---

## Full Workflow

This repo is the second step in a two-pipeline workflow:

```
Step 1: bovine-igg-pipeline
  Raw nanopore reads → demultiplexing → filtering → consensus sequences
  https://github.com/fruggles11/bovine-igg-pipeline
  Output: results/4_consensus_sequences/  (and optionally 5_majority_consensus/)

Step 2: bovine-repertoire-analysis (this repo)
  Consensus sequences → V(D)J annotation → diversity analysis + ultra-long CDR3H3 detection
```

Run `bovine-igg-pipeline` with `--skip_annotation true` and do all V(D)J/CDR3 calling here instead. Its own optional annotation step (`--skip_annotation false`) calls IgBLAST with native AIRR output (`-outfmt 19`) and has no junction-inference fallback — without proper auxiliary data for the custom bovine J germlines, `junction`/`junction_aa`/`productive` come back blank for essentially every sequence. This pipeline's `INFER_JUNCTION` step (see [CDR3/Junction Calling](#cdrjunction-calling)) avoids that failure mode.

### Recommended workflow for ultra-long CDR3H3 detection

If `main.nf` is running correctly:

```bash
# Run main repertoire pipeline (Step 2)
nextflow run . \
    --germline_db './bovine_germline' \
    --input_dir '/path/to/results/4_consensus_sequences' \
    --results ./repertoire_results

# Then extract ultra-long CDR3H3 sequences from the AIRR output
python3 find_ultralong_cdrh3.py \
    --airr_dir repertoire_results/filtered \
    --output_dir ultralong_cdrh3
```

If `main.nf` is unavailable, use `ultralong_cdrh3_filter.nf` as a fallback (see [below](#ultra-long-cdrh3-filter-ultralong_cdrh3_filternf)):

```bash
# Run the standalone ultra-long CDR3H3 pipeline instead
nextflow run fruggles11/bovine-repertoire-analysis -main-script ultralong_cdrh3_filter.nf \
    --fasta_dir /path/to/results/5_majority_consensus \
    --germline_db './bovine_germline'

# Then extract ultra-long CDR3H3 sequences
python3 find_ultralong_cdrh3.py \
    --airr_dir ultralong_results/filtered \
    --fasta_dir /path/to/results/5_majority_consensus \
    --germline bovine_germline/Bos_taurus_IgHJ_gaps.fasta \
    --output_dir ultralong_cdrh3
```

---

## Pipelines

| Pipeline | Script | Purpose |
|----------|--------|---------|
| Repertoire analysis | `main.nf` | Full V(D)J annotation, clonotype assignment, diversity analysis |
| Ultra-long CDR3H3 filter | `ultralong_cdrh3_filter.nf` | Fallback — runs IgBLAST on majority consensus FASTAs only, without full diversity analysis |

## Features

- **V(D)J annotation** using IgBLAST
- **CDR3/junction fallback calling** — backfills CDR3 calls IgBLAST's heuristic misses via J-gene anchor scanning (see [CDR3/Junction Calling](#cdrjunction-calling))
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
# Bovine germline sequences are bundled in the bovine_germline/ directory.
# Clone the repo and point --germline_db at it:
git clone https://github.com/fruggles11/bovine-repertoire-analysis
cd bovine-repertoire-analysis

nextflow run . \
    --germline_db './bovine_germline' \
    --input_dir '/path/to/results/4_consensus_sequences' \
    --results ./repertoire_results
```

Or run directly from GitHub (germlines must be downloaded separately — see [Germline Database](#germline-database-required)):

```bash
nextflow run fruggles11/bovine-repertoire-analysis \
    --germline_db '/path/to/germlines' \
    --input_dir '/path/to/results/4_consensus_sequences' \
    --results ./repertoire_results
```

The pipeline will automatically find all `*.fasta` files recursively within the input directory.

## Input

The pipeline accepts FASTA files containing antibody sequences. These can be:
- `*_unique.fasta` files from the bovine-igg-pipeline
- `*_majority.fasta` files from the bovine-igg-pipeline `5_majority_consensus/` step
- `*_nogroup_unique.fasta` files
- Any FASTA with antibody V(D)J sequences

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_dir` | none | Directory to recursively search for `*.fasta` files (recommended) |
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
├── igblast/                # Raw IgBLAST annotations
├── airr/                   # AIRR-formatted sequence data (direct from MakeDb.py)
├── airr_junction_filled/   # Same, with missing CDR3/junction calls backfilled
│                           # (see CDR3/Junction Calling below)
├── filtered/                # Productive sequences only
├── clones/                  # Clone assignments
├── reports/
│   └── vdj_summary.tsv     # One row per sequence: barcode, chain, v/d/j_call,
│                           # productive, v_identity, junction_aa_length, junction_source
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

The pipeline requires bovine immunoglobulin germline sequences. A bundled copy is included in `bovine_germline/` — point `--germline_db` at this directory when running locally.

If you need to re-download from IMGT, **you must download IMGT-gapped sequences.** Ungapped sequences will not work correctly — CDR3/junction regions will not be identified and all downstream analysis will fail.

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

---

## CDR3/Junction Calling

IgBLAST/Change-O's `MakeDb.py` sometimes leaves `junction`/`junction_aa` blank for a sequence even though it successfully assigned V/D/J genes — the pipeline's `INFER_JUNCTION` step (`main.nf`, using `bin/infer_missing_junction.py`) backfills these gaps.

For heavy chain sequences (`IGHV` in `v_call`) with no `junction_aa`, it scans the sequence for the conserved J-region Trp codon (matched against the IGHJ germline) and walks back to the nearest upstream Cys codon — the same J-gene anchor method used by `find_ultralong_cdrh3.py`. Unlike that script, this fallback has no minimum length floor, since the goal is calling CDR3 length generally rather than just flagging ultra-long sequences. It only fills gaps — it never overwrites a junction IgBLAST already called.

Every row gets a `junction_source` column:
- `igblast` — IgBLAST/MakeDb.py called the junction directly
- `anchor` — filled in by the J-gene anchor fallback
- blank — neither method could confidently call it (e.g. ambiguous bases, non-standard 3′ ends, or a light-chain sequence contaminating a heavy-chain sample — the fallback is heavy-chain only, since light chain CDR3s are short and reliably called by IgBLAST already)

In practice, most gaps turn out to be light-chain sequences misassigned into heavy-chain barcode files (see [Automatic Chain Filtering](#automatic-chain-filtering)) rather than genuine IgBLAST failures — check `v_call` before assuming a blank junction is a detection bug.

---

## Ultra-Long CDR3H3 Filter (`ultralong_cdrh3_filter.nf`)

> **Note:** This pipeline is a lightweight fallback. If `main.nf` is available, use that instead — it produces full diversity analysis and AIRR output, and `find_ultralong_cdrh3.py` can be run on its `filtered/` output.

Bovine heavy chains uniquely produce CDR3H3 regions up to 70+ amino acids long — far exceeding the ~15 aa typical in humans. This pipeline runs IgBLAST directly on the majority consensus FASTAs from [bovine-igg-pipeline](https://github.com/fruggles11/bovine-igg-pipeline) to identify these ultra-long sequences.

Only sequences with a confirmed IGHV gene call are reported — IgBLAST occasionally assigns light chain germlines (IGLV) to heavy chain sequences, so this filter prevents false positives.

### Quick Start

```bash
nextflow run fruggles11/bovine-repertoire-analysis -main-script ultralong_cdrh3_filter.nf \
    --fasta_dir /path/to/results/5_majority_consensus \
    --germline_db './bovine_germline'
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

### Notes on CDR3H3 detection

IgBLAST reports `Total identifiable CDR3 = 0` for some ultra-long sequences because the CDR3 length exceeds its internal heuristic window. In practice this means sequences with very long CDR3s (>55 aa) are reliably detected, while sequences in the 40–50 aa range may occasionally be missed or have slightly underestimated lengths. Running with `--min_cdr3_aa 40` captures the full ultra-long population with minimal false positives when the IGHV filter is applied.

### Utilities

**`find_ultralong_cdrh3.py`** — standalone script that combines IgBLAST-based and J gene anchor CDR3H3 detection. Run this after either `main.nf` or `ultralong_cdrh3_filter.nf` to extract the final ultra-long CDR3H3 hit list. IgBLAST results take precedence (more accurate lengths); the J gene anchor method runs as a fallback on barcodes IgBLAST did not detect.

```bash
# After main.nf
python3 find_ultralong_cdrh3.py \
    --airr_dir repertoire_results/filtered \
    --output_dir ultralong_cdrh3

# After ultralong_cdrh3_filter.nf (adds J anchor fallback)
python3 find_ultralong_cdrh3.py \
    --airr_dir  ultralong_results/filtered \
    --fasta_dir results/5_majority_consensus \
    --germline  bovine_germline/Bos_taurus_IgHJ_gaps.fasta \
    --output_dir ultralong_cdrh3
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--airr_dir` | required | Directory with `*_heavy_productive.tsv` AIRR files |
| `--fasta_dir` | optional | Directory with `*_heavy_majority.fasta` files — enables J gene anchor fallback |
| `--germline` | optional | Bovine IGHJ germline FASTA — required when `--fasta_dir` is provided |
| `--output_dir` | `ultralong_cdrh3` | Output directory |
| `--min_cdr3_aa` | `40` | Minimum CDR3H3 length in amino acids |
| `--min_j_identity` | `0.75` | Minimum J gene anchor match identity |

Output files in `--output_dir`:
- `ultralong_cdrh3.tsv` — annotation table with CDR3H3 length, junction sequence, V/J/D calls, and `source` column (`igblast` or `anchor`)
- `ultralong_cdrh3.fasta` — sequences sorted by CDR3H3 length descending; headers include `v_call`, `cdrh3_aa`, and `source`

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
