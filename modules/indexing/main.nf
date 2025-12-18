process buildMinimapIndexGenome {
    cpus 8
    memory 64.GB
    time '12h'

    input:
    path(ref_genome)

    output:
    path(index_mmi)

    script:
    index_mmi = "${ref_genome}.mmi"
    """
    minimap2 -x splice:hq -d ${index_mmi} ${ref_genome}
    """
}

process buildMinimapIndexTranscriptome {
    cpus 8
    memory 64.GB
    time '12h'

    input:
    path(ref)

    output:
    path(index_mmi)

    script:
    index_mmi = "${ref}.mmi"
    """
    minimap2 -x lr:hq -d ${index_mmi} ${ref}
    """
}
