def getBarcodeUmiRegex(chemistry) {
    // Barcode/UMI lengths per 10X Chromium chemistry, used to build the nailpolish
    // header regex. Must stay consistent with getChemistryPatterns in
    // modules/barcode_detection/main.nf (barcode/umi '?' counts).
    def lengths = [
        '3v1'  : [bc: 14, umi: 10],
        '3v2'  : [bc: 16, umi: 10],
        '3v3'  : [bc: 16, umi: 12],
        '3v3.1': [bc: 16, umi: 12],
        '3v4'  : [bc: 16, umi: 12],
        '5v1'  : [bc: 16, umi: 10],
        '5v2'  : [bc: 16, umi: 10],
        '5v3'  : [bc: 16, umi: 12],
    ]
    def dims = lengths[chemistry]
    if (dims == null) {
        error "Unsupported chemistry version: ${chemistry}. Supported versions: ${lengths.keySet().join(', ')}"
    }
    return "^(?<CB>[ATCGNX]{${dims.bc}})_(?<UB>[ATCGNX]{${dims.umi}})"
}

process nailpolishDedup {
    label 'large'
    tag "$sample_id"
    publishDir "${params.output_dir}/dedup_fastq/",
        mode: 'copy',
        pattern: '*_dedup.fastq.gz',
        enabled: params.publish_dedup_fastq ?: false
    publishDir "${params.output_dir}/qc/nailpolish/",
        mode: 'copy',
        pattern: '*.summary.html',
        enabled: params.publish_qc_dedup_logs ?: false

    input:
    tuple val(sample_id), path(tagged_fastq)

    output:
    tuple val(sample_id), path(dedup_fastq), emit: fastq
    path(summary_html),                      emit: summary

    script:
    // Molecular (CB+UMI) consensus deduplication at the FASTQ level. Reads are grouped
    // by the barcode/UMI encoded in the flexiplex read header; false duplicates are
    // split by the clustering step, and each true group is consensus-called (error
    // corrected) into a single read. CB:Z/UB:Z tags are carried in the FASTQ comment
    // and propagated into the BAM by minimap2 -y downstream.
    def barcode_regex = getBarcodeUmiRegex(params.chemistry)
    reads_fastq  = "${sample_id}_tagged.fastq"
    dedup_fastq  = "${sample_id}_dedup.fastq.gz"
    summary_html = "${sample_id}.summary.html"
    """
    gunzip -c ${tagged_fastq} > ${reads_fastq}

    nailpolish index ${reads_fastq} \\
        --barcode-regex "${barcode_regex}" \\
        --skip-unmatched

    nailpolish summary ${reads_fastq} -o ${summary_html}

    nailpolish consensus ${reads_fastq} \\
        -t ${task.cpus} \\
        -o ${sample_id}_dedup.fastq

    pigz -p ${task.cpus} ${sample_id}_dedup.fastq
    """
}

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
        --per-cell \\
        --cell-tag=CB \\
        --extract-umi-method=tag \\
        --umi-tag=UB \\
        -L ${dedup_log} \\
        -I ${input_bam} \\
        -S ${dedup_bam}
    """
}
