process flexiplexGetBarcodeCandidates {
    publishDir "${params.output_dir}/flexiplex/",
        mode: 'copy',
        pattern: '*.txt',
        enabled: params.publish.flexiplex_candidates
    tag "${sample_id}"
    cpus 8
    memory 16.GB
    time '4h'

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path(output_barcodes)

    script:
    output_barcodes = "${sample_id}_barcodes_counts.txt"
    """
    gunzip -c ${fastq_file} | flexiplex -p ${task.cpus} -x "CTACACGACGCTCTTCCGATCT" -b "????????????????" -u "????????????" -x "TTTCTTATATGGG" -f 0 -n ${sample_id}
    """
}

process mergeFlexiplexBarcodes {
    publishDir "${params.output_dir}/flexiplex/merged/",
        mode: 'copy',
        enabled: params.publish.flexiplex_merged
    cpus 4
    memory 8.GB
    time '2h'

    input:
    path(barcode_files)
    path(canonical_bc_list)

    output:
    path(merged_bc_file)

    script:
    merged_bc_file = "flexiplex_merged_barcodes.txt"
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    bc_list <- read_lines("${canonical_bc_list}")
    files <- fs::dir_ls(".", glob = "*barcodes_counts.txt")
    data_raw <- map(files, ~read_tsv(., col_names = c("barcode", "count")))

    merge_counts <- function(x, y) {
        full_join(
            rename(x, count_x = count),
            rename(y, count_y = count)
        ) %>%
            replace_na(list(count_x = 0, count_y = 0)) %>%
            mutate(count = count_x + count_y, .keep = "unused")
    }

    data <- reduce(data_raw, merge_counts) %>%
        filter(barcode %in% bc_list) %>%
        filter(count > 500)

    data %>%
        select(barcode) %>%
        write_tsv('${merged_bc_file}', col_names = FALSE)

    """
}

process flexiplexTagFastq {
    publishDir "${params.output_dir}/flexiplex/barcoded_fastq/",
        mode: 'copy',
        pattern: '*.fastq.gz',
        enabled: params.publish.barcoded_fastq
    publishDir "logs/flexiplex",
        mode: 'copy',
        pattern: "*.log",
        enabled: params.publish.flexiplex_logs
    tag "${sample_id}"
    cpus 8
    memory 32.GB
    time '48h'

    input:
    tuple val(sample_id), path(fastq_file), path(barcode_file)

    output:
    tuple val(sample_id), path(output_fastq), emit: 'fastq'
    path(output_log), emit: 'log'

    script:
    output_fastq = "${sample_id}_tagged.fastq.gz"
    output_log = "${sample_id}_flexiplex.log"
    """
    gunzip -c ${fastq_file} | flexiplex -p ${task.cpus} -x "CTACACGACGCTCTTCCGATCT" -b "????????????????" -u "????????????" -x "TTTCTTATATGGG" -f 2 -e 1 -n ${sample_id} -k ${barcode_file} | pigz > ${output_fastq}
    grep -vE '^INFO:|^WARNING:|million reads processed|cite us' .command.log > ${output_log}
    """
}
