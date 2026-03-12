process runBcl2Fastq {
    // /stornext/System/data/software/rhel/9/base/bioinf/bcl2fastq/2.20.0.422/bin/bcl2fastq
    label 'medium'
    publishDir "${params.output_dir}/bcl2fastq/",
        mode: 'copy',
        enabled: params.publish_bcl2fastq

    input:
    path(bcl_folder)

    output:
    path(output_fastq)

    script:
    output_fastq = "hashtag_fastq"
    """
    bcl2fastq \
        --runfolder-dir ${bcl_folder} \
        --output-dir ${output_fastq} \
        --no-lane-splitting
    """
}
