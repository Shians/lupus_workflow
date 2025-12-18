process bamToFastq {
    publishDir "${params.output_dir}/fastq/",
        mode: 'copy',
        pattern: '*.fastq.gz',
        enabled: params.publish.fastq
    cpus 4
    memory 16.GB
    time '12h'
    array 100

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path(fastq_file)

    script:
    fastq_file = "${sample_id}.fastq.gz"
    """
    samtools fastq -F 0x900 -@ ${task.cpus} ${bam_file} | pigz > ${fastq_file}
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

process retagFastq {
    publishDir "${params.output_dir}/retagged_fastq/",
        mode: 'copy',
        pattern: '*.fastq.gz',
        enabled: params.publish.retagged_fastq
    tag "${sample_id}"
    cpus 8
    memory 16.GB
    time '48h'

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path(output_fastq)

    script:
    output_fastq = "${sample_id}_retagged.fastq.gz"
    """
    gunzip -c ${fastq_file} | awk '
        BEGIN {OFS=""}
        {
            line = NR % 4
            if (line == 1) {
                # Example: @ATTCGGGGTGTGTCGA_CGCAACACTGCT#867acfcb-f5d4-4c70-a950-aa1ce68bdc1e_-1of1
                header = substr(\$0, 2)  # remove leading @
                n = split(header, parts, "#")
                if (n == 2) {
                    id = parts[1]        # barcode_umi
                    readname = parts[2]  # readname
                    split(id, id_parts, "_")
                    barcode = id_parts[1]
                    umi = id_parts[2]
                    print "@" id "#" readname "\tCB:Z:" barcode "\tUB:Z:" umi
                } else {
                    print \$0  # fallback
                }
            } else {
                print \$0
            }
        }
    ' | pigz > ${output_fastq}
    """
}
