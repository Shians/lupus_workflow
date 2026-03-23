process umitoolsDedup {
    label 'medium'
    tag "$sample_id"
    publishDir "${params.output_dir}/bam_dedup/",
        mode: 'copy',
        pattern: '*.{bam,bai}',
        enabled: params.publish_bam_dedup ?: false

    input:
    tuple val(sample_id), path(input_bam), path(input_bai)

    output:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai), emit: bam

    script:
    dedup_bam = "${sample_id}_dedup.bam"
    dedup_bai = "${dedup_bam}.bai"
    """
    umi_tools dedup \\
        --per-contig \\
        --per-cell \\
        --cell-tag=CB \\
        --extract-umi-method=tag \\
        --umi-tag=UB \\
        -I ${input_bam} \\
        -S ${dedup_bam}
    samtools index ${dedup_bam}
    """
}
