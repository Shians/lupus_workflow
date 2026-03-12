process runCellSNPGenotype {
    label 'large'
    publishDir "${params.output_dir}/cell_snp/",
        mode: 'copy',
        enabled: params.publish_cell_snp
    tag "CellSNP-genotype ${suffix}"

    input:
    tuple path(bam_paths), path(path_indices), path(barcode_path), val(suffix), path(snp_annotation)

    output:
    tuple path(output_path), val(suffix)

    script:
    output_path = "cellsnp_" + suffix
    bam_paths = bam_paths.join(",")
    """
    echo ${bam_paths}
    mkdir -p cellsnp
    cellsnp-lite -s ${bam_paths} \
        -b ${barcode_path} \
        -O ${output_path} \
        -R ${snp_annotation} \
        -p ${task.cpus} \
        --minMAF 0.1 \
        --minCOUNT 10 \
        --gzip --genotype
    """
}

process runVireoDemultiplex {
    label 'large'
    publishDir "${params.output_dir}/vireo/",
        mode: 'copy',
        enabled: params.publish_vireo
    tag "Vireo-demultiplex ${suffix}"

    input:
    tuple path(cellsnp_path), val(suffix), val(n_donors)

    output:
    path output_path

    script:
    output_path = "vireo_" + suffix
    """
    mkdir -p vireo
    vireo \
        -c ${cellsnp_path} \
        -o ${output_path} \
        -N ${n_donors} \
        -t GT \
        -p ${task.cpus}
    """
}
