process indexBam {
    label 'medium'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path(input_bam), path(input_bai), emit: bam_bai

    script:
    input_bai = "${input_bam}.bai"
    """
    samtools index ${input_bam}
    """
}

process sortIndexBam {
    label 'medium'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path(sorted_bam), path(sorted_bai)

    script:
    sorted_bam = "${sample_id}_coord_sorted.bam"
    sorted_bai = "${sorted_bam}.bai"
    """
    samtools sort -@ ${task.cpus} -o ${sorted_bam} ${input_bam}
    samtools index ${sorted_bam}
    """
}

process umitoolsDedup {
    label 'large'
    tag "$sample_id"
    publishDir "${params.output_dir}/bam_dedup/",
        mode: 'copy',
        pattern: '*.bam',
        enabled: params.publish_bam_dedup ?: false
    publishDir "${params.output_dir}/logs/umitools/",
        mode: 'copy',
        pattern: '*.log',
        enabled: params.publish_qc_dedup_logs ?: false

    input:
    tuple val(sample_id), path(input_bam), path(input_bai)

    output:
    tuple val(sample_id), path(dedup_bam), emit: bam
    path(dedup_log), emit: log

    script:
    dedup_bam = "${sample_id}_dedup.bam"
    dedup_log = "${sample_id}_dedup.log"
    """
    umi_tools dedup \\
        --per-gene \\
        --per-contig \\
        --per-cell \\
        --cell-tag=CB \\
        --extract-umi-method=tag \\
        --umi-tag=UB \\
        -L ${dedup_log} \\
        -I ${input_bam} \\
        -S ${dedup_bam}
    """
}
