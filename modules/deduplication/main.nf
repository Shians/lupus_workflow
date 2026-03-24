process umitoolsDedup {
    label 'medium'
    tag "$sample_id"
    publishDir "${params.output_dir}/bam_dedup/",
        mode: 'copy',
        pattern: '*.bam',
        enabled: params.publish_bam_dedup ?: false

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path(dedup_bam), emit: bam

    script:
    dedup_bam = "${sample_id}_dedup.bam"
    """
    samtools index ${input_bam}
    umi_tools dedup \\
        --per-gene \\
        --per-contig \\
        --per-cell \\
        --cell-tag=CB \\
        --extract-umi-method=tag \\
        --umi-tag=UB \\
        -I ${input_bam} \\
        -S ${dedup_bam}
    """
}
