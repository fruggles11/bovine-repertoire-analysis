#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    Bovine IgG Repertoire Analysis Pipeline
========================================================================================
    Analyzes antibody repertoire diversity using the Immcantation framework

    Tools used:
    - IgBLAST: V(D)J gene assignment
    - Change-O: Clonal assignment and AIRR formatting
    - Alakazam: Diversity analysis and visualization

    Input: FASTA files from bovine-igg-pipeline (unique sequences)
    Output: Clonotype diversity metrics, CDR3 analysis, gene usage plots
========================================================================================
*/

// --------------------------------------------------------------- //
// WORKFLOW
// --------------------------------------------------------------- //

workflow {

    // Input channel from FASTA files
    ch_fasta = Channel
        .fromPath( params.fasta_input )
        .map { fasta ->
            def sample_id = fasta.baseName.replaceAll(/_unique$/, '').replaceAll(/_nogroup$/, '')
            tuple( sample_id, fasta )
        }

    // Germline database is required for bovine (not available via fetch_imgtdb.sh)
    if ( !params.germline_db ) {
        error """
        =====================================================================
        ERROR: Bovine germline database is required
        =====================================================================
        The Immcantation fetch_imgtdb.sh script does not support bovine.
        You must provide germline FASTA files manually.

        Download from IMGT/GENE-DB (https://www.imgt.org/genedb/):
          1. Select Species: "Bos taurus"
          2. Select Gene type: IGHV, IGHD, IGHJ (and IGKV, IGKJ, IGLV, IGLJ if needed)
          3. Download nucleotide sequences in FASTA format
          4. Save files to a directory (e.g., ./germlines/)

        Then run with:
          nextflow run fruggles11/bovine-repertoire-analysis \\
              --germline_db './germlines/*.fasta' \\
              --fasta_input 'analysis_input/*_unique.fasta'
        =====================================================================
        """.stripIndent()
    }
    ch_germline = Channel.fromPath( params.germline_db )

    // Build IgBLAST database
    BUILD_IGBLAST_DB( ch_germline.collect() )

    // Run IgBLAST annotation
    IGBLAST_ANNOTATION(
        ch_fasta,
        BUILD_IGBLAST_DB.out.db,
        BUILD_IGBLAST_DB.out.aux
    )

    // Convert IgBLAST output to AIRR format
    MAKEDB(
        IGBLAST_ANNOTATION.out,
        ch_germline.collect()
    )

    // Filter for productive sequences
    FILTER_PRODUCTIVE(
        MAKEDB.out
    )

    // Define clones
    DEFINE_CLONES(
        FILTER_PRODUCTIVE.out
    )

    // Collect all samples for combined analysis
    ch_all_clones = DEFINE_CLONES.out.collect()

    // Run diversity analysis
    DIVERSITY_ANALYSIS(
        ch_all_clones
    )

    // Generate summary report
    GENERATE_REPORT(
        DIVERSITY_ANALYSIS.out.plots,
        DIVERSITY_ANALYSIS.out.stats
    )

}

// --------------------------------------------------------------- //
// PARAMETERS
// --------------------------------------------------------------- //

params.fasta_input = "${launchDir}/*_unique.fasta"
params.germline_db = null  // Path to pre-built germline FASTA files
params.species = "bovine"
params.chain = "IG"  // IG for immunoglobulin
params.results = "${launchDir}/results"

// Clone definition parameters
params.clone_threshold = 0.15  // Junction hamming distance threshold
params.clone_method = "hierarchical"  // clustering method

// Analysis parameters
params.min_seq_count = 1  // Minimum sequences per clone
params.subsample = null  // Subsample size for diversity (null = no subsampling)
params.skip_productive_filter = true  // Productivity detection not working reliably for bovine - skipping filter

// Barcode/chain assignment for filtering cross-contamination
// Specify which barcodes contain only light or heavy chains
// Any other chain type found in these barcodes will be filtered out as contamination
params.light_chain_barcodes = "barcode88"  // Comma-separated list of barcodes that are light chain only
params.heavy_chain_barcodes = "barcode96"  // Comma-separated list of barcodes that are heavy chain only


// --------------------------------------------------------------- //
// PROCESSES
// --------------------------------------------------------------- //

// NOTE: FETCH_GERMLINES is not used for bovine - germlines must be provided manually
// This process is kept for potential future use with other species
process FETCH_GERMLINES {

    publishDir "${params.results}/germlines", mode: 'copy'

    output:
    path "*.fasta"

    script:
    """
    # fetch_imgtdb.sh only supports: human, mouse, rat, rabbit, rhesus_monkey
    # Bovine germlines must be downloaded manually from IMGT
    echo "ERROR: Automatic germline fetch not supported for bovine"
    echo "Please download from https://www.imgt.org/genedb/ and use --germline_db"
    exit 1
    """
}

process BUILD_IGBLAST_DB {

    publishDir "${params.results}/igblast_db", mode: 'copy'

    input:
    path germlines

    output:
    path "database/*", emit: db
    path "internal_data/*", emit: aux

    script:
    """
    mkdir -p database internal_data

    # Debug: show what germline files we received
    echo "Germline files received:" >&2
    ls -la >&2

    # Separate V, D, J genes based on filename patterns
    # Handle filenames with spaces by using find and proper quoting
    # Use -iname for case-insensitive matching (handles IGHV, IgHV, ighv, etc.)

    # V genes (IGHV, IGKV, IGLV)
    find . -maxdepth 1 -iname "*IGHV*" -o -iname "*IGKV*" -o -iname "*IGLV*" -o -iname "*IgHV*" -o -iname "*IgKV*" -o -iname "*IgLV*" | while read f; do
        cat "\$f" >> database/bovine_V.fasta 2>/dev/null
    done

    # D genes (IGHD only - light chains don't have D)
    find . -maxdepth 1 -iname "*IGHD*" -o -iname "*IgHD*" | while read f; do
        cat "\$f" >> database/bovine_D.fasta 2>/dev/null
    done

    # J genes (IGHJ, IGKJ, IGLJ)
    find . -maxdepth 1 -iname "*IGHJ*" -o -iname "*IGKJ*" -o -iname "*IGLJ*" -o -iname "*IgHJ*" -o -iname "*IgKJ*" -o -iname "*IgLJ*" | while read f; do
        cat "\$f" >> database/bovine_J.fasta 2>/dev/null
    done

    # Ensure files exist (even if empty)
    touch database/bovine_V.fasta database/bovine_D.fasta database/bovine_J.fasta

    # Debug: show what we created
    echo "Database files created:" >&2
    ls -la database/ >&2
    echo "V file lines: \$(wc -l < database/bovine_V.fasta 2>/dev/null || echo 0)" >&2
    echo "D file lines: \$(wc -l < database/bovine_D.fasta 2>/dev/null || echo 0)" >&2
    echo "J file lines: \$(wc -l < database/bovine_J.fasta 2>/dev/null || echo 0)" >&2

    # Build BLAST databases only if files have content
    cd database
    for f in bovine_*.fasta; do
        if [[ -s "\$f" ]]; then
            echo "Building BLAST database for \$f" >&2
            makeblastdb -parse_seqids -dbtype nucl -in "\$f"
        else
            echo "Warning: \$f is empty, skipping" >&2
        fi
    done
    cd ..

    # Create IgBLAST auxiliary data
    # Copy internal data from IgBLAST installation if IGDATA is set
    if [[ -n "\${IGDATA:-}" ]] && [[ -d "\${IGDATA}/internal_data" ]]; then
        cp -r "\${IGDATA}/internal_data/"* internal_data/ 2>/dev/null || true
    fi

    # Create organism-specific files if not present
    mkdir -p internal_data/bovine
    # Create basic aux file for bovine (required for IgBLAST)
    cat > internal_data/bovine/bovine_gl.aux << 'AUXFILE'
# Bovine germline auxiliary data
# Frame information for bovine IG genes
AUXFILE
    """
}

process IGBLAST_ANNOTATION {

    tag "${sample_id}"
    publishDir "${params.results}/igblast", mode: 'copy'

    cpus 4

    input:
    tuple val(sample_id), path(fasta)
    path db
    path aux

    output:
    tuple val(sample_id), path("${sample_id}_igblast.fmt7"), path(fasta)

    script:
    """
    # Set up IgBLAST directory structure
    mkdir -p igblast_data/database

    # Copy our bovine database files (FASTA and BLAST index files)
    # The db variable contains all files from BUILD_IGBLAST_DB output
    for f in ${db}; do
        cp "\$f" igblast_data/database/
    done

    # Debug: list what we copied
    echo "Database files copied:" >&2
    ls -la igblast_data/database/ >&2

    # Copy entire internal_data from container's IGDATA (includes human databases IgBLAST needs)
    if [[ -n "\${IGDATA:-}" ]] && [[ -d "\${IGDATA}/internal_data" ]]; then
        cp -r "\${IGDATA}/internal_data" igblast_data/
    else
        # Fallback: try common container paths
        for path in /usr/local/share/igblast/internal_data /usr/share/igblast/internal_data; do
            if [[ -d "\$path" ]]; then
                cp -r "\$path" igblast_data/
                break
            fi
        done
    fi

    # Set IGDATA to our custom directory
    export IGDATA="\$(pwd)/igblast_data"

    # Run IgBLAST directly (AssignGenes.py doesn't support bovine)
    # Use -organism human for internal_data validation (required by IgBLAST)
    # The actual annotation uses our bovine databases via -germline_db_* flags
    # Note: BLAST db names include .fasta because makeblastdb was run on .fasta files
    igblastn \
        -query ${fasta} \
        -out ${sample_id}_igblast.fmt7 \
        -num_threads ${task.cpus} \
        -ig_seqtype Ig \
        -organism human \
        -germline_db_V "\${IGDATA}/database/bovine_V.fasta" \
        -germline_db_D "\${IGDATA}/database/bovine_D.fasta" \
        -germline_db_J "\${IGDATA}/database/bovine_J.fasta" \
        -outfmt "7 std qseq sseq btop" \
        -domain_system imgt
    """
}

process MAKEDB {

    tag "${sample_id}"
    publishDir "${params.results}/airr", mode: 'copy'

    input:
    tuple val(sample_id), path(igblast_out), path(fasta)
    path germlines

    output:
    tuple val(sample_id), path("${sample_id}_db-pass.tsv")

    script:
    """
    # Combine germline references using find to handle spaces in filenames
    # Use -iname for case-insensitive matching (handles IGHV, IgHV, ighv, etc.)
    find . -maxdepth 1 -iname "*IGHV*" -o -iname "*IGKV*" -o -iname "*IGLV*" -o -iname "*IgHV*" -o -iname "*IgKV*" -o -iname "*IgLV*" | while read f; do
        cat "\$f" >> combined_V.fasta 2>/dev/null
    done
    find . -maxdepth 1 -iname "*IGHD*" -o -iname "*IgHD*" | while read f; do
        cat "\$f" >> combined_D.fasta 2>/dev/null
    done
    find . -maxdepth 1 -iname "*IGHJ*" -o -iname "*IGKJ*" -o -iname "*IGLJ*" -o -iname "*IgHJ*" -o -iname "*IgKJ*" -o -iname "*IgLJ*" | while read f; do
        cat "\$f" >> combined_J.fasta 2>/dev/null
    done

    # Ensure files exist
    touch combined_V.fasta combined_D.fasta combined_J.fasta

    # Debug: show germline file sizes
    echo "Germline files for MakeDb:" >&2
    ls -la combined_*.fasta >&2

    # Convert IgBLAST output to AIRR format
    # Use --partial to allow sequences without complete V(D)J
    # Use --infer-junction to extract junction/CDR3 from alignments (needed for bovine)
    MakeDb.py igblast \
        -i ${igblast_out} \
        -s ${fasta} \
        -r combined_V.fasta combined_D.fasta combined_J.fasta \
        --extended \
        --partial \
        --infer-junction \
        --format airr \
        -o ${sample_id}_db.tsv

    # Rename output to expected filename
    # MakeDb.py with -o sample_db.tsv outputs to sample_db.tsv (not sample_db_db-pass.tsv)
    if [[ -f "${sample_id}_db.tsv" ]] && [[ ! -f "${sample_id}_db-pass.tsv" ]]; then
        mv ${sample_id}_db.tsv ${sample_id}_db-pass.tsv
    elif [[ -f "${sample_id}_db_db-pass.tsv" ]]; then
        mv ${sample_id}_db_db-pass.tsv ${sample_id}_db-pass.tsv
    fi

    # If no pass file was created, create empty one to not fail pipeline
    if [[ ! -f "${sample_id}_db-pass.tsv" ]]; then
        echo "Warning: No sequences passed MakeDb, creating empty file" >&2
        echo -e "sequence_id\\tv_call\\td_call\\tj_call\\tsequence\\tproductive" > ${sample_id}_db-pass.tsv
    fi
    """
}

process FILTER_PRODUCTIVE {

    tag "${sample_id}"
    publishDir "${params.results}/filtered", mode: 'copy'

    input:
    tuple val(sample_id), path(airr_tsv)

    output:
    tuple val(sample_id), path("${sample_id}_productive.tsv")

    script:
    if (params.skip_productive_filter)
        """
        # Productivity filtering skipped (--skip_productive_filter true)
        # Use --skip_productive_filter false with IMGT-gapped germlines
        echo "Skipping productivity filter (germlines not IMGT-gapped)" >&2
        cp ${airr_tsv} ${sample_id}_productive.tsv
        """
    else
        """
        # Filter for productive sequences (in-frame, no stop codons)
        # Requires IMGT-gapped germline sequences for accurate results
        line_count=\$(wc -l < ${airr_tsv})

        if [[ \$line_count -le 1 ]]; then
            echo "Warning: Input file is empty or header-only, copying as-is" >&2
            cp ${airr_tsv} ${sample_id}_productive.tsv
        else
            ParseDb.py select \
                -d ${airr_tsv} \
                -f productive \
                -u T TRUE True true \
                -o ${sample_id}_productive.tsv || cp ${airr_tsv} ${sample_id}_productive.tsv
        fi

        # Ensure output exists
        if [[ ! -f "${sample_id}_productive.tsv" ]]; then
            cp ${airr_tsv} ${sample_id}_productive.tsv
        fi
        """
}

process DEFINE_CLONES {

    tag "${sample_id}"
    publishDir "${params.results}/clones", mode: 'copy'

    input:
    tuple val(sample_id), path(airr_tsv)

    output:
    path "${sample_id}_clones.tsv"

    script:
    """
    # Check if file has junction field and has data
    line_count=\$(wc -l < ${airr_tsv})
    has_junction=\$(head -1 ${airr_tsv} | grep -c "junction" || echo "0")

    if [[ \$line_count -le 1 ]] || [[ \$has_junction -eq 0 ]]; then
        echo "Warning: Input file lacks junction field or is empty, copying as-is" >&2
        cp ${airr_tsv} ${sample_id}_clones.tsv
    else
        # Define clones based on V gene, J gene, and junction similarity
        DefineClones.py -d ${airr_tsv} \
            --act set \
            --model ham \
            --norm len \
            --dist ${params.clone_threshold} \
            --format airr \
            -o ${sample_id}_clones.tsv || cp ${airr_tsv} ${sample_id}_clones.tsv
    fi

    # Ensure output exists
    if [[ ! -f "${sample_id}_clones.tsv" ]]; then
        cp ${airr_tsv} ${sample_id}_clones.tsv
    fi
    """
}

process DIVERSITY_ANALYSIS {

    publishDir "${params.results}/diversity", mode: 'copy'

    input:
    path clone_files

    output:
    path "plots/*", emit: plots
    path "stats/*", emit: stats

    script:
    """
    mkdir -p plots stats

    # Run R script for diversity analysis
    Rscript - <<'RSCRIPT'

    library(alakazam)
    library(shazam)
    library(ggplot2)
    library(dplyr)
    library(airr)

    # Read all clone files
    files <- list.files(pattern = "_clones.tsv\$", full.names = TRUE)

    if (length(files) == 0) {
        stop("No clone files found")
    }

    # Combine all data
    db_list <- lapply(files, function(f) {
        df <- read_rearrangement(f)
        # Get base sample name, removing _clones.tsv suffix
        sample_name <- gsub("_clones.tsv", "", basename(f))
        # Remove _0, _1, _2, etc. and _nogroup suffixes to combine related samples
        # e.g., barcode96_heavy_0, barcode96_heavy_1 -> barcode96_heavy
        sample_name <- gsub("_[0-9]+\$", "", sample_name)
        sample_name <- gsub("_nogroup\$", "", sample_name)
        df\$sample_id <- sample_name
        return(df)
    })
    db <- bind_rows(db_list)

    # Filter out cross-contamination based on expected barcode/chain combinations
    # Light chain barcodes should not have heavy chain data (and vice versa)
    light_barcodes <- strsplit("${params.light_chain_barcodes}", ",")[[1]]
    light_barcodes <- trimws(light_barcodes)
    heavy_barcodes <- strsplit("${params.heavy_chain_barcodes}", ",")[[1]]
    heavy_barcodes <- trimws(heavy_barcodes)

    # Build contamination list: light barcodes + _heavy, heavy barcodes + _light
    contamination <- c(
        paste0(light_barcodes, "_heavy"),
        paste0(heavy_barcodes, "_light")
    )
    contamination <- contamination[contamination != "_heavy" & contamination != "_light"]  # Remove empty entries

    if (length(contamination) > 0) {
        db <- db %>% filter(!sample_id %in% contamination)
        message("Filtered out contamination samples: ", paste(contamination, collapse = ", "))
    }
    message("Remaining samples: ", paste(unique(db\$sample_id), collapse = ", "))

    # Check if clone_id column exists
    has_clone_id <- "clone_id" %in% colnames(db)
    has_junction_aa <- "junction_aa" %in% colnames(db)

    if (!has_clone_id) {
        message("Warning: clone_id column not found in data. Clone-based analyses will be skipped.")
    }

    # Basic stats
    stats <- db %>%
        group_by(sample_id) %>%
        summarise(
            total_sequences = n(),
            unique_clones = if (has_clone_id) n_distinct(clone_id, na.rm = TRUE) else NA_integer_,
            productive = sum(productive == TRUE | productive == "T", na.rm = TRUE),
            mean_cdr3_length = if (has_junction_aa) mean(nchar(as.character(junction_aa)), na.rm = TRUE) else NA_real_,
            median_cdr3_length = if (has_junction_aa) median(nchar(as.character(junction_aa)), na.rm = TRUE) else NA_real_
        )
    write.csv(stats, "stats/basic_stats.csv", row.names = FALSE)

    # CDR3 length distribution (only if junction_aa exists)
    if (has_junction_aa) {
        db\$cdr3_length <- nchar(as.character(db\$junction_aa))

        p1 <- ggplot(db, aes(x = cdr3_length, fill = sample_id)) +
            geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
            labs(title = "CDR3 Length Distribution",
                 x = "CDR3 Length (amino acids)",
                 y = "Count") +
            theme_minimal() +
            theme(legend.position = "bottom")
        ggsave("plots/cdr3_length_distribution.pdf", p1, width = 10, height = 6)
        ggsave("plots/cdr3_length_distribution.png", p1, width = 10, height = 6, dpi = 150)
    } else {
        message("Skipping CDR3 length distribution due to missing junction_aa column")
    }

    # V gene usage
    v_usage <- db %>%
        filter(!is.na(v_call)) %>%
        mutate(v_gene = gsub("\\\\*.*", "", v_call)) %>%
        group_by(sample_id, v_gene) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(sample_id) %>%
        mutate(freq = count / sum(count))

    write.csv(v_usage, "stats/v_gene_usage.csv", row.names = FALSE)

    p2 <- ggplot(v_usage, aes(x = reorder(v_gene, -freq), y = freq, fill = sample_id)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "V Gene Usage",
             x = "V Gene",
             y = "Frequency") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
              legend.position = "bottom")
    ggsave("plots/v_gene_usage.pdf", p2, width = 14, height = 6)
    ggsave("plots/v_gene_usage.png", p2, width = 14, height = 6, dpi = 150)

    # J gene usage
    j_usage <- db %>%
        filter(!is.na(j_call)) %>%
        mutate(j_gene = gsub("\\\\*.*", "", j_call)) %>%
        group_by(sample_id, j_gene) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(sample_id) %>%
        mutate(freq = count / sum(count))

    write.csv(j_usage, "stats/j_gene_usage.csv", row.names = FALSE)

    p3 <- ggplot(j_usage, aes(x = reorder(j_gene, -freq), y = freq, fill = sample_id)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "J Gene Usage",
             x = "J Gene",
             y = "Frequency") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")
    ggsave("plots/j_gene_usage.pdf", p3, width = 10, height = 6)
    ggsave("plots/j_gene_usage.png", p3, width = 10, height = 6, dpi = 150)

    # Clone size distribution (only if clone_id exists)
    if (has_clone_id) {
        clone_sizes <- db %>%
            filter(!is.na(clone_id)) %>%
            group_by(sample_id, clone_id) %>%
            summarise(clone_size = n(), .groups = "drop")

        write.csv(clone_sizes, "stats/clone_sizes.csv", row.names = FALSE)

        if (nrow(clone_sizes) > 0) {
            p4 <- ggplot(clone_sizes, aes(x = clone_size, fill = sample_id)) +
                geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
                scale_x_log10() +
                labs(title = "Clone Size Distribution",
                     x = "Clone Size (log10)",
                     y = "Count") +
                theme_minimal() +
                theme(legend.position = "bottom")
            ggsave("plots/clone_size_distribution.pdf", p4, width = 10, height = 6)
            ggsave("plots/clone_size_distribution.png", p4, width = 10, height = 6, dpi = 150)
        }

        # Diversity curves using Alakazam
        if (nrow(db) > 0 && any(!is.na(db\$clone_id))) {

            # Calculate diversity
            div_curve <- tryCatch({
                alphaDiversity(db, group = "sample_id", clone = "clone_id",
                              min_q = 0, max_q = 4, step_q = 0.1,
                              ci = 0.95, nboot = 100)
            }, error = function(e) {
                message("Could not calculate diversity curve: ", e\$message)
                NULL
            })

            if (!is.null(div_curve)) {
                p5 <- plot(div_curve, legend_title = "Sample") +
                    labs(title = "Repertoire Diversity (Hill Numbers)") +
                    theme_minimal()
                ggsave("plots/diversity_curve.pdf", p5, width = 10, height = 6)
                ggsave("plots/diversity_curve.png", p5, width = 10, height = 6, dpi = 150)

                # Save diversity values
                write.csv(div_curve@diversity, "stats/diversity_values.csv", row.names = FALSE)
            }

            # Rarefaction curve
            rarefaction <- tryCatch({
                estimateAbundance(db, group = "sample_id", clone = "clone_id",
                                ci = 0.95, nboot = 100)
            }, error = function(e) {
                message("Could not calculate rarefaction: ", e\$message)
                NULL
            })

            if (!is.null(rarefaction)) {
                p6 <- plot(rarefaction, legend_title = "Sample") +
                    labs(title = "Clonal Abundance Rarefaction") +
                    theme_minimal()
                ggsave("plots/rarefaction_curve.pdf", p6, width = 10, height = 6)
                ggsave("plots/rarefaction_curve.png", p6, width = 10, height = 6, dpi = 150)
            }
        }

        # Summary diversity metrics
        diversity_summary <- db %>%
            filter(!is.na(clone_id)) %>%
            group_by(sample_id) %>%
            summarise(
                total_sequences = n(),
                unique_clones = n_distinct(clone_id),
                simpson_index = 1 - sum((table(clone_id)/n())^2),
                shannon_index = -sum((table(clone_id)/n()) * log(table(clone_id)/n())),
                chao1 = n_distinct(clone_id) + (sum(table(clone_id) == 1)^2) / (2 * max(1, sum(table(clone_id) == 2))),
                .groups = "drop"
            )

        write.csv(diversity_summary, "stats/diversity_summary.csv", row.names = FALSE)
    } else {
        message("Skipping clone-based analyses due to missing clone_id column")
        # Create empty placeholder files
        write.csv(data.frame(), "stats/clone_sizes.csv", row.names = FALSE)
        write.csv(data.frame(), "stats/diversity_summary.csv", row.names = FALSE)
    }

    message("Diversity analysis complete!")

    RSCRIPT
    """
}

process GENERATE_REPORT {

    publishDir "${params.results}", mode: 'copy'

    input:
    path plots
    path stats

    output:
    path "repertoire_report.html"

    script:
    """
    # Generate HTML report
    cat > repertoire_report.html << 'HTML'
<!DOCTYPE html>
<html>
<head>
    <title>Bovine IgG Repertoire Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2c3e50; }
        h2 { color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        .stats-table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        .stats-table th, .stats-table td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        .stats-table th { background-color: #3498db; color: white; }
        .stats-table tr:nth-child(even) { background-color: #f2f2f2; }
        .plot-container { margin: 20px 0; text-align: center; }
        .plot-container img { max-width: 100%; border: 1px solid #ddd; }
        .section { margin: 30px 0; }
    </style>
</head>
<body>
    <h1>Bovine IgG Repertoire Analysis Report</h1>
    <p>Generated: \$(date)</p>

    <div class="section">
        <h2>Summary Statistics</h2>
        <p>See stats/ directory for detailed CSV files.</p>
    </div>

    <div class="section">
        <h2>CDR3 Length Distribution</h2>
        <div class="plot-container">
            <img src="diversity/plots/cdr3_length_distribution.png" alt="CDR3 Length Distribution">
        </div>
    </div>

    <div class="section">
        <h2>V Gene Usage</h2>
        <div class="plot-container">
            <img src="diversity/plots/v_gene_usage.png" alt="V Gene Usage">
        </div>
    </div>

    <div class="section">
        <h2>J Gene Usage</h2>
        <div class="plot-container">
            <img src="diversity/plots/j_gene_usage.png" alt="J Gene Usage">
        </div>
    </div>

    <div class="section">
        <h2>Clone Size Distribution</h2>
        <div class="plot-container">
            <img src="diversity/plots/clone_size_distribution.png" alt="Clone Size Distribution">
        </div>
    </div>

    <div class="section">
        <h2>Diversity Analysis</h2>
        <div class="plot-container">
            <img src="diversity/plots/diversity_curve.png" alt="Diversity Curve">
        </div>
        <div class="plot-container">
            <img src="diversity/plots/rarefaction_curve.png" alt="Rarefaction Curve">
        </div>
    </div>

    <div class="section">
        <h2>Output Files</h2>
        <ul>
            <li><strong>airr/</strong> - AIRR-formatted sequence annotations</li>
            <li><strong>clones/</strong> - Clone assignments</li>
            <li><strong>diversity/stats/</strong> - CSV files with diversity metrics</li>
            <li><strong>diversity/plots/</strong> - PDF and PNG plots</li>
        </ul>
    </div>

</body>
</html>
HTML
    """
}
