process catTranscriptAlignedBams {
    tag "$sample_id"
    cpus 8
    memory 32.GB
    time '12h'

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
    publishDir "${params.output_dir}/splice_aligned/",
        mode: 'copy',
        pattern: '*.{bam,bai}',
        enabled: params.publish.splice_aligned
    tag "$sample_id"
    cpus 8
    memory 32.GB
    time '12h'

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

process combineMergedSplicedBams {
    cpus 8
    memory 32.GB
    time '12h'

    input:
    path(bam_files, stageAs: "*.bam")

    output:
    tuple path(merged_bam), path(merged_bam_index)

    script:
    merged_bam = "merged.bam"
    merged_bam_index = "${merged_bam}.bai"
    """
    samtools merge -@ ${task.cpus} -o ${merged_bam} ${bam_files}
    samtools index ${merged_bam}
    """
}

process sortBamByName {
    tag "$sample"
    cpus 8
    memory "16.GB"

    input:
    tuple val(sample), path(input_bam)

    output:
    tuple val(sample), path(output_bam)

    script:
    output_bam = "${sample}_sorted.bam"
    """
    samtools sort -@ ${task.cpus} -n -o ${output_bam} ${input_bam}
    """
}
