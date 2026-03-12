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

    // Set barcode and UMI lengths based on chemistry
    def chemistry = params.chemistry
    def barcode_pattern = ""
    def umi_pattern = ""
    def suffix_pattern = ""

    if (chemistry == '3v1') {
        barcode_pattern = "??????????????"     // 14 bp
        umi_pattern = "??????????"             // 10 bp
        suffix_pattern = "TTTTTTTT"
    } else if (chemistry == '3v2') {
        barcode_pattern = "????????????????"   // 16 bp
        umi_pattern = "??????????"             // 10 bp
        suffix_pattern = "TTTTTTTT"
    } else if (chemistry in ['3v3', '3v3.1', '3v4']) {
        barcode_pattern = "????????????????"   // 16 bp
        umi_pattern = "????????????"           // 12 bp
        suffix_pattern = "TTTTTTTT"
    } else if (chemistry in ['5v1', '5v2']) {
        barcode_pattern = "????????????????"   // 16 bp
        umi_pattern = "??????????"             // 10 bp
        suffix_pattern = "TTTCTTATAT"
    } else if (chemistry == '5v3') {
        barcode_pattern = "????????????????"   // 16 bp
        umi_pattern = "????????????"           // 12 bp
        suffix_pattern = "TTTCTTATAT"
    } else {
        error "Unsupported chemistry version: ${chemistry}. Supported versions: 3v1, 3v2, 3v3, 3v3.1, 3v4, 5v1, 5v2, 5v3"
    }

    """
    # full left-pattern is CTACACGACGCTCTTCCGATCT, but we can allow for some mismatches in the first 20 bp to capture more reads
    gunzip -c ${fastq_file} | flexiplex -p ${task.cpus} -x "GACGCTCTTCCGATCT" -b "${barcode_pattern}" -u "${umi_pattern}" -x "${suffix_pattern}" -f 0 -n ${sample_id}
    """
}

process mergeFlexiplexBarcodes {
    label 'small'
    publishDir "${params.output_dir}/flexiplex_merged_barcodes/",
        mode: 'copy',
        enabled: params.publish_flexiplex_merged
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(barcode_files, stageAs: '*_barcodes_counts.txt'), path(canonical_bc_list)

    output:
    tuple val(sample_id), path(merged_bc_file)

    script:
    merged_bc_file = "${sample_id}_flexiplex_merged_barcodes.txt"
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

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
    """
}

process flexiplexTagFastq {
    label 'large'
    publishDir "${params.output_dir}/flexiplex_barcoded_fastq/",
        mode: 'copy',
        pattern: '*.fastq.gz',
        enabled: params.publish_retagged_fastq
    publishDir "logs/flexiplex",
        mode: 'copy',
        pattern: "*.log",
        enabled: params.publish_flexiplex_logs
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_file), path(barcode_file)

    output:
    tuple val(sample_id), path(output_fastq), emit: 'fastq'
    path(output_log), emit: 'log'

    script:
    output_fastq = "${sample_id}_tagged.fastq.gz"
    output_log = "${sample_id}_flexiplex.log"

    // Set barcode and UMI lengths based on chemistry
    def chemistry = params.chemistry
    def barcode_pattern = ""
    def umi_pattern = ""
    def suffix_pattern = ""

    if (chemistry == '3v1') {
        barcode_pattern = "??????????????"     // 14 bp
        umi_pattern = "??????????"             // 10 bp
        suffix_pattern = "TTTTTTTT"
    } else if (chemistry == '3v2') {
        barcode_pattern = "????????????????"   // 16 bp
        umi_pattern = "??????????"             // 10 bp
        suffix_pattern = "TTTTTTTT"
    } else if (chemistry in ['3v3', '3v3.1', '3v4']) {
        barcode_pattern = "????????????????"   // 16 bp
        umi_pattern = "????????????"           // 12 bp
        suffix_pattern = "TTTTTTTT"
    } else if (chemistry in ['5v1', '5v2']) {
        barcode_pattern = "????????????????"   // 16 bp
        umi_pattern = "??????????"             // 10 bp
        suffix_pattern = "TTTCTTATAT"
    } else if (chemistry == '5v3') {
        barcode_pattern = "????????????????"   // 16 bp
        umi_pattern = "????????????"           // 12 bp
        suffix_pattern = "TTTCTTATAT"
    } else {
        error "Unsupported chemistry version: ${chemistry}. Supported versions: 3v1, 3v2, 3v3, 3v3.1, 3v4, 5v1, 5v2, 5v3"
    }

    """
    gunzip -c ${fastq_file} | flexiplex -p ${task.cpus} -x "GACGCTCTTCCGATCT" -b "${barcode_pattern}" -u "${umi_pattern}" -x "${suffix_pattern}" -f 2 -e 1 -n ${sample_id} -k ${barcode_file} | pigz > ${output_fastq}
    grep -vE '^INFO:|^WARNING:|million reads processed|cite us' .command.log > ${output_log}
    """
}
