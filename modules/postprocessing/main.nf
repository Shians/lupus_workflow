process catTranscriptAlignedBams {
    label 'medium'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam_files, stageAs: "*.bam")

    output:
    tuple val(sample_id), path(merged_bam)

    script:
    merged_bam = "${sample_id}.bam"
    """
    samtools cat -o ${merged_bam} ${bam_files}
    """
}

process mergeSpliceAlignedBams {
    label 'medium'
    publishDir "${params.output_dir}/splice_aligned/",
        mode: 'copy',
        pattern: '*.{bam,bai}',
        enabled: params.publish_splice_aligned
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam_files, stageAs: "*.bam")

    output:
    tuple val(sample_id), path(merged_bam), path(merged_bam_index)

    script:
    merged_bam = "${sample_id}.bam"
    merged_bam_index = "${merged_bam}.bai"
    """
    samtools merge -@ ${task.cpus} -o ${merged_bam} ${bam_files}
    samtools index ${merged_bam}
    """
}

process sortBamByCellBarcode {
    label 'medium'
    tag "$sample"

    input:
    tuple val(sample), path(input_bam)

    output:
    tuple val(sample), path(output_bam)

    script:
    output_bam = "${sample}_cb_sorted.bam"
    """
    samtools sort -@ ${task.cpus} -t CB -o ${output_bam} ${input_bam}
    """
}
