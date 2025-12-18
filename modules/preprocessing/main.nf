process bamToFastq {
    publishDir "${params.output_dir}/fastq/",
        mode: 'copy',
        pattern: '*.fastq.gz',
        enabled: params.publish.fastq
    tag "${sample_id}-${bam_file.simpleName}"
    cpus 4
    memory '16.GB'
    time '12h'
    array 100

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path(fastq_file)

    script:
    fastq_file = "${bam_file.simpleName}.fastq.gz"
    """
    samtools fastq -F 0x900 -@ ${task.cpus} ${bam_file} | pigz -p ${task.cpus} > ${fastq_file}
    """
}

process splitFastqChunks {
    tag "$sample"
    cpus 8
    memory "16.GB"

    input:
    tuple val(sample), path(fastq), val(alignment_bam_parts)

    output:
    tuple val(sample), path("chunks/${sample}.part_*.fastq.gz")

    script:
    """
    mkdir -p chunks
    seqkit split2 \
        -f \
        -O chunks \
        --threads ${task.cpus} \
        --by-part ${alignment_bam_parts} \
        --by-part-prefix ${sample}.part_ \
        --extension ".gz" \
        --compress-level 1 \
        $fastq
    """
}

