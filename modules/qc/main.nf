process craminoStats {
    label 'tiny'
    tag "${sample_id}-${bam_file.simpleName}"
    publishDir "${params.output_dir}/qc/cramino/",
        mode: 'copy',
        enabled: params.publish_qc_cramino ?: false

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path(cramino_out)

    script:
    cramino_out = "${sample_id}_${bam_file.simpleName}_cramino.json"
    """
    cramino --ubam --format json ${bam_file} > ${cramino_out}
    """
}

process samtoolsFlagstat {
    label 'tiny'
    tag "${sample_id}"
    publishDir "${params.output_dir}/qc/flagstat/",
        mode: 'copy',
        enabled: params.publish_qc_flagstat ?: false

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path(flagstat_file)

    script:
    flagstat_file = "${sample_id}_flagstat.txt"
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${flagstat_file}
    """
}

process mosdepthCoverage {
    label 'small'
    tag "${sample_id}"
    publishDir "${params.output_dir}/qc/mosdepth/",
        mode: 'copy',
        enabled: params.publish_qc_mosdepth ?: false

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id),
        path("${sample_id}.mosdepth.summary.txt"),
        path("${sample_id}.mosdepth.global.dist.txt")

    script:
    """
    mosdepth --fast-mode -t ${task.cpus} ${sample_id} ${bam}
    """
}

process countReads {
    label 'tiny'
    tag "${sample_id} (${stage})"

    input:
    tuple val(sample_id), path(bam), val(stage), val(count_mode)

    output:
    tuple val(sample_id), path(count_file)

    script:
    count_file = "${sample_id}_${stage}_count.txt"
    def mode_flags = [
        ALL            : '',
        PRIMARY        : '-F 0x900',
        PRIMARY_MAPPED : '-F 0x904'
    ]
    def flags = mode_flags[count_mode]
    """
    count=\$(samtools view -c ${flags} ${bam})
    echo -e "${sample_id}\t${stage}\t\$count" > ${count_file}
    """
}

process countFastqReads {
    label 'tiny'
    tag "${sample_id} (${stage})"

    input:
    tuple val(sample_id), path(fastq_file), val(stage)

    output:
    tuple val(sample_id), path(count_file)

    script:
    // FASTQ counterpart to countReads, for stages that sit between tools rather
    // than between BAMs. Every read in a FASTQ is a primary read by definition,
    // so there is no count_mode to select here.
    count_file = "${sample_id}_${stage}_count.txt"
    """
    lines=\$(pigz -dc -p ${task.cpus} ${fastq_file} | wc -l)
    echo -e "${sample_id}\t${stage}\t\$(( lines / 4 ))" > ${count_file}
    """
}

process readCountSummary {
    label 'tiny'
    tag "${sample_id}"
    publishDir "${params.output_dir}/qc/read_tracking/",
        mode: 'copy',
        enabled: params.publish_qc_read_tracking ?: false

    input:
    tuple val(sample_id), path(count_files, stageAs: '*_count.txt')

    output:
    path(summary_tsv)

    script:
    summary_tsv = "${sample_id}_read_tracking.tsv"
    """
    echo -e "sample_id\tstage\treads" > ${summary_tsv}
    cat *_count.txt >> ${summary_tsv}
    """
}

process aggregateReadTracking {
    label 'tiny'
    publishDir "${params.output_dir}/qc/read_tracking/",
        mode: 'copy',
        enabled: params.publish_qc_read_tracking ?: false

    input:
    tuple path(tracking_files, stageAs: 'tracking/*'), val(stage_order)

    output:
    path(combined_tsv), emit: summary
    path(multiqc_tsv),  emit: multiqc

    script:
    // stage_order is a comma-separated "stage:branch" list built in the entry
    // workflow, since which stages exist depends on params.dedup_method. It
    // fixes row order (the per-sample files arrive in nondeterministic task
    // completion order) and tells the fraction maths which stages sit on the
    // shared trunk vs. on the genome/transcriptome forks.
    combined_tsv = "read_tracking_all_samples.tsv"
    multiqc_tsv  = "read_tracking_mqc.tsv"
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    stage_spec <- str_split_1("${stage_order}", ",")
    stage_meta <- tibble(
        stage       = str_split_i(stage_spec, ":", 1),
        branch      = str_split_i(stage_spec, ":", 2),
        stage_index = seq_along(stage_spec)
    )

    # Counts arrive one file per input BAM, not one per sample (bamToFastq and
    # flexiplexTagFastq both fan out over the BAMs listed for a sample), so the
    # per-stage totals have to be summed before any fraction is meaningful.
    counts <- list.files("tracking", pattern = "[.]tsv\$", full.names = TRUE) %>%
        map(read_tsv, show_col_types = FALSE) %>%
        bind_rows() %>%
        summarise(reads = sum(reads), .by = c(sample_id, stage))

    tracked <- counts %>%
        inner_join(stage_meta, by = "stage") %>%
        arrange(sample_id, stage_index)

    add_fractions <- function(sample_rows) {
        input_reads <- first(pull(sample_rows, reads))
        # Both forks descend from the last shared stage, so that is the
        # denominator for the first row of each fork.
        trunk_reads <- sample_rows %>%
            filter(branch == "common") %>%
            slice_tail(n = 1) %>%
            pull(reads)

        sample_rows %>%
            group_by(branch) %>%
            mutate(prev_reads = lag(reads)) %>%
            ungroup() %>%
            mutate(
                prev_reads = if_else(
                    is.na(prev_reads) & branch != "common", trunk_reads, prev_reads
                ),
                pct_of_input    = round(100 * reads / input_reads, 2),
                pct_of_previous = round(100 * reads / prev_reads, 2)
            )
    }

    summary_tbl <- tracked %>%
        group_split(sample_id) %>%
        map(add_fractions) %>%
        bind_rows() %>%
        select(sample_id, branch, stage, reads, pct_of_input, pct_of_previous)

    write_tsv(summary_tbl, "${combined_tsv}")

    # MultiQC custom content. A table rather than a bargraph: the per-stage
    # counts are nested subsets of one another, so stacking them would imply an
    # additive breakdown that does not exist.
    mqc_header <- c(
        "# id: 'read_tracking'",
        "# section_name: 'Read Tracking'",
        "# description: 'Primary reads surviving each pipeline stage, summed over all input BAMs per sample.'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'read_tracking_table'",
        "#     title: 'Primary reads by pipeline stage'"
    )

    mqc_table <- summary_tbl %>%
        select(sample_id, stage, reads) %>%
        pivot_wider(names_from = stage, values_from = reads)

    write_lines(mqc_header, "${multiqc_tsv}")
    write_tsv(mqc_table, "${multiqc_tsv}", append = TRUE, col_names = TRUE)
    """
}

process multiQC {
    label 'small'
    publishDir "${params.output_dir}/qc/",
        mode: 'copy',
        enabled: params.publish_qc_multiqc ?: false

    input:
    path(qc_files)

    output:
    path("multiqc_report.html")
    path("multiqc_data/")

    script:
    """
    multiqc .
    """
}
