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
    // nailpolish can read gzipped FASTQ, but doing so incurs a large performance hit
    // (every index/summary/consensus pass re-inflates the file), so the tagged
    // .fastq.gz is decompressed once up front and all commands run on the plain FASTQ.
    def barcode_regex = getBarcodeUmiRegex(params.chemistry)
    tagged_plain = "${sample_id}_tagged.fastq"
    dedup_fastq  = "${sample_id}_dedup.fastq.gz"
    summary_html = "${sample_id}.summary.html"
    """
    pigz -dc -p ${task.cpus} ${tagged_fastq} > ${tagged_plain}

    nailpolish index ${tagged_plain} \\
        --barcode-regex "${barcode_regex}" \\
        --skip-unmatched

    nailpolish summary ${tagged_plain} -o ${summary_html}

    nailpolish consensus ${tagged_plain} \\
        -t ${task.cpus} \\
        -o ${sample_id}_dedup.fastq

    pigz -p ${task.cpus} ${sample_id}_dedup.fastq
    """
}
