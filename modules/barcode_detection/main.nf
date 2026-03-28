def getChemistryPatterns(chemistry) {
    if (chemistry == '3v1') {
        return [prefix: "GACGCTCTTCCGATCT", barcode: "??????????????", umi: "??????????", suffix: "TTTTTTTT"]
    } else if (chemistry == '3v2') {
        return [prefix: "GACGCTCTTCCGATCT", barcode: "????????????????", umi: "??????????", suffix: "TTTTTTTT"]
    } else if (['3v3', '3v3.1', '3v4'].contains(chemistry)) {
        return [prefix: "GACGCTCTTCCGATCT", barcode: "????????????????", umi: "????????????", suffix: "TTTTTTTT"]
    } else if (['5v1', '5v2'].contains(chemistry)) {
        return [prefix: "GACGCTCTTCCGATCT", barcode: "????????????????", umi: "??????????", suffix: "TTTCTTATAT"]
    } else if (chemistry == '5v3') {
        return [prefix: "GACGCTCTTCCGATCT", barcode: "????????????????", umi: "????????????", suffix: "TTTCTTATAT"]
    } else {
        error "Unsupported chemistry version: ${chemistry}. Supported versions: 3v1, 3v2, 3v3, 3v3.1, 3v4, 5v1, 5v2, 5v3"
    }
}

process flexiplexGetBarcodeCandidates {
    label 'medium'
    publishDir "${params.output_dir}/flexiplex/",
        mode: 'copy',
        pattern: '*.txt',
        enabled: params.publish_flexiplex_candidates
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path(output_barcodes)

    script:
    output_barcodes = "${sample_id}_barcodes_counts.txt"

    def chemistryPatterns = getChemistryPatterns(params.chemistry)

    """
    gunzip -c ${fastq_file} \\
        | flexiplex \\
            -p ${task.cpus} \\
            -x "${chemistryPatterns.prefix}" \\
            -b "${chemistryPatterns.barcode}" \\
            -u "${chemistryPatterns.umi}" \\
            -x "${chemistryPatterns.suffix}" \\
            -f 0 \\
            -n ${sample_id}
    """
}

process mergeFlexiplexBarcodes {
    label 'small'
    publishDir "${params.output_dir}/flexiplex_merged_barcodes/",
        mode: 'copy',
        pattern: '*.txt',
        enabled: params.publish_flexiplex_merged
    publishDir "${params.output_dir}/qc/knee_plots/",
        mode: 'copy',
        pattern: '*.png',
        enabled: params.publish_qc_knee_plots ?: false
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(barcode_files, stageAs: '*_barcodes_counts.txt'), path(canonical_bc_list)

    output:
    tuple val(sample_id), path(merged_bc_file), emit: barcodes
    path(knee_plot), emit: knee_plot

    script:
    merged_bc_file = "${sample_id}_flexiplex_merged_barcodes.txt"
    knee_plot = "${sample_id}_barcode_knee_plot.png"
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(ggplot2)

    bc_list <- read_lines("${canonical_bc_list}")

    # Read all barcode count files and combine them
    barcode_files <- list.files(pattern = "*_barcodes_counts.txt")

    data <- map(barcode_files, ~read_tsv(.x, col_names = c("barcode", "count"), show_col_types = FALSE)) %>%
        bind_rows() %>%
        summarize(count = sum(count), .by = barcode) %>%
        arrange(desc(count))

    data %>%
        filter(barcode %in% bc_list) %>%
        filter(count > ${params.min_barcode_count}) %>%
        select(barcode) %>%
        write_tsv('${merged_bc_file}', col_names = FALSE)

    plot_data <- data %>%
        filter(barcode %in% bc_list) %>%
        mutate(rank = row_number())

    n_total <- nrow(plot_data)
    n_passed <- sum(plot_data\$count > ${params.min_barcode_count})
    pct <- ifelse(n_total > 0, round(n_passed / n_total * 100, 1), 0)
    subtitle <- glue::glue(
        "{n_passed} / {n_total} barcodes passed ({pct}%),  min count = {min_count}",
        min_count = ${params.min_barcode_count}
    )

    ggplot(plot_data, aes(x = rank, y = count)) +
        geom_line() +
        geom_vline(
            xintercept = n_passed,
            linetype = "dashed", colour = "red"
        ) +
        scale_x_log10() +
        scale_y_log10() +
        labs(
            title = "${sample_id} — Barcode Knee Plot",
            subtitle = subtitle,
            x = "Barcode rank",
            y = "Read count"
        ) +
        theme_bw()

    ggsave("${knee_plot}", width = 7, height = 5)
    """
}

process flexiplexTagFastq {
    label 'large'
    publishDir "${params.output_dir}/flexiplex_barcoded_fastq/",
        mode: 'copy',
        pattern: '*.fastq.gz|*.txt.gz',
        enabled: params.publish_retagged_fastq
    publishDir "${params.output_dir}/logs/flexiplex",
        mode: 'copy',
        pattern: "*.log",
        enabled: params.publish_flexiplex_logs
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_file), path(barcode_file)

    output:
    tuple val(sample_id), path(output_fastq),          emit: fastq
    tuple val(sample_id), path(output_reads_barcodes), emit: reads_barcodes
    path(output_log),                                  emit: log

    script:
    output_fastq = "${sample_id}_${fastq_file.simpleName}_tagged.fastq.gz"
    output_reads_barcodes = "${sample_id}_reads_barcodes.txt.gz"
    output_log = "${sample_id}_${fastq_file.simpleName}_flexiplex.log"

    def chemistryPatterns = getChemistryPatterns(params.chemistry)

    """
    gunzip -c ${fastq_file} \\
        | flexiplex \\
            -p ${task.cpus} \\
            -x "${chemistryPatterns.prefix}" \\
            -b "${chemistryPatterns.barcode}" \\
            -u "${chemistryPatterns.umi}" \\
            -x "${chemistryPatterns.suffix}" \\
            -f 2 \\
            -e 1 \\
            -n ${sample_id} \\
            -k ${barcode_file} \\
        | pigz > ${output_fastq}
    pigz ${sample_id}_reads_barcodes.txt
    grep -vE '^INFO:|^WARNING:|million reads processed|cite us' .command.log > ${output_log}
    """
}
