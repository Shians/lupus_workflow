process runOarfish {
    label 'large'
    publishDir "${params.output_dir}/oarfish/",
        mode: "copy",
        enabled: params.publish_oarfish
    tag "$sample"
    errorStrategy { task.exitStatus in [137, 139] ? 'retry' : 'finish' }
    maxRetries 3

    input:
    tuple val(sample), path(input_bam)

    output:
    tuple val(sample), path(output_dir)

    script:
    output_dir = "oarfish_output/${sample}"
    """
    oarfish \
        -j ${task.cpus} \
        --single-cell \
        --model-coverage \
        --filter-group no-filters \
        --output ${output_dir}/sample \
        --alignments ${input_bam}
    """
}
