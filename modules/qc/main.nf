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
    mosdepth --fast-mode --no-abbrev -t ${task.cpus} ${sample_id} ${bam}
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
        PRIMARY_MAPPED : '-F 0x904'
    ]
    def flags = mode_flags[count_mode]
    """
    count=\$(samtools view -c ${flags} ${bam})
    echo -e "${sample_id}\t${stage}\t\$count" > ${count_file}
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
