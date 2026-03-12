process alignMinimap2Spliced {
    label 'large'
    tag "$sample"
    array 1000

    input:
    tuple val(sample), path(fastq), path(ref)

    output:
    tuple val(sample), path(out_bam)

    script:
    out_bam = fastq.baseName + ".bam"

    """
    minimap2 -t ${task.cpus} -ax splice:hq -ub -y ${ref} ${fastq} | \
        samtools view -b | \
        samtools sort > ${out_bam}
    samtools index ${out_bam}
    """
}

process alignMinimap2TranscriptomeUnsorted {
    label 'large'
    tag "$sample"
    array 1000

    input:
    tuple val(sample), path(fastq), path(ref)

    output:
    tuple val(sample), path(out_bam)

    script:
    out_bam = fastq.baseName + ".bam"

    """
    minimap2 -t ${task.cpus} -a -y --rev-only ${ref} ${fastq} | \
        samtools view -b | \
        samtools sort > ${out_bam}
    samtools index ${out_bam}
    """
}
