// cellsnp-lite site filters scaled to the size of the donor pool.
//
// minMAF: cellsnp-lite pools REF/ALT counts across every cell in the library,
//   so for D equally represented donors a SNP where exactly one donor is
//   heterozygous -- the most abundant informative class, and the one that most
//   cleanly fingerprints an individual donor -- has pooled MAF 0.5/D. Filtering
//   at 0.5/D would sit exactly on that boundary, where uneven donor
//   representation or sampling noise drops the site, so the threshold is halved
//   to 0.25/D. That tolerates a donor contributing as little as half its
//   balanced share of cells. Capped at 0.1 (the conventional small-pool value;
//   a looser filter buys nothing when donors are few) and floored at 0.01
//   (below that, residual ONT error competes with true signal even on the
//   curated -T list).
//
// minCOUNT: derived from minMAF rather than from the donor count directly. At
//   a site with N molecules the observable MAFs are 0, 1/N, 2/N, ..., so a
//   minMAF below 1/N passes every site carrying even one minor-allele molecule
//   and the filter silently stops discriminating. Requiring N >= 1/minMAF makes
//   the threshold mean "at least one minor-allele molecule", the weakest form
//   that still does work. Floored at the cellsnp-lite/Vireo standard of 20,
//   which is a data-quality floor and does not depend on D.
//
//   Note this does NOT scale because Vireo needs per-donor depth at each site.
//   It does not: Vireo pools evidence across the whole matrix, and assignment
//   aggregates many sites per cell rather than resolving each site across all
//   donors. Larger pools want more sites, not deeper ones. If a large-D run
//   collapses the retained site count on ONT depth, raise the minMAF floor
//   below and let minCOUNT follow it down -- do not chase depth this data
//   does not have.
//
// Shared by the monolithic and sharded processes so the two paths cannot drift
// apart. Both filters are computed per-SNP from that SNP's own pileup, so they
// are shard-invariant and the sharded output stays comparable to a single run.
def cellsnpDonorFilters(n_donors) {
    def donors = n_donors as Integer
    if (donors < 1) {
        error "n_donors must be >= 1 to scale cellSNP filters, got: ${n_donors}"
    }
    def min_maf = Math.min(0.1d, Math.max(0.01d, 0.25d / donors))
    def min_count = Math.max(20, Math.ceil(1.0d / min_maf) as Integer)
    return [minMAF: String.format('%.4f', min_maf), minCOUNT: min_count]
}

// Normalise the target SNP list to coordinate order once, before it fans out
// to either cellSNP path.
//
// cellsnp-lite emits sites in target-file order, so the order of this file is
// the order of every downstream cellSNP directory. Sorting here means the
// monolithic and sharded paths are comparable row-for-row rather than only as
// sets, and it is the precondition for ever merging sharded chunks back into
// coordinate order: a round-robin chunk of a sorted VCF is itself sorted, so
// the merge can be an ordered k-way merge instead of a blind concatenation.
//
// Sorting only. The VCF is taken as given -- a missing ##fileformat line,
// absent ##contig headers, or duplicate records are the caller's to fix, and
// bcftools will complain about them far more clearly than we could.
process sortSNPAnnotation {
    label 'small'
    tag "Sort SNP annotation"

    input:
    path snp_annotation

    output:
    path output_path, emit: vcf
    path contig_order_path, emit: contig_order

    script:
    output_path = "sorted_snp_annotation.vcf.gz"
    contig_order_path = "contig_order.txt"
    // Leave headroom below the tier limit: -m caps bcftools' in-core buffer,
    // not the process, and exceeding it spills to -T rather than failing.
    def sort_memory = Math.max(1, (task.memory.toGiga() as Integer) - 2)
    """
    bcftools sort \
        -m ${sort_memory}G \
        -T . \
        -Oz \
        -o ${output_path} \
        ${snp_annotation}

    # The contig ordering bcftools actually produced, which mergeCellSNP needs
    # to rebuild that same order from the chunks. Taken from the records rather
    # than the ##contig header lines, which need not be present or complete.
    # Sorting makes each contig contiguous, so uniq collapses it to first
    # appearance.
    bcftools query -f '%CHROM\\n' ${output_path} | uniq > ${contig_order_path}
    """
}

// Cheap depth sieve ahead of cellSNP genotyping.
//
// mosdepth --by only needs alignment coordinates, not per-base quality or a
// genotype-likelihood model, so it is dramatically cheaper per site than
// cellsnp-lite's own per-cell pileup: this touches each sample's BAM once,
// where cellsnp-lite would touch it once per cell barcode. Sites with no
// coverage in a sample contribute nothing to that sample's genotyping and
// would be dropped by cellsnp-lite's own minCOUNT filter anyway (see
// cellsnpDonorFilters above) -- this just does that filtering before the
// expensive per-cell pass instead of after it.
process computeTargetDepth {
    label 'medium'
    tag "Target-depth ${sample_id}"

    input:
    tuple path(bam), path(bai), val(sample_id), path(snp_annotation)

    output:
    path output_path

    script:
    output_path = "depth_${sample_id}.bed"
    """
    bcftools query -f '%CHROM\\t%POS\\n' ${snp_annotation} \
        | awk -v OFS='\t' '{print \$1, \$2-1, \$2}' > candidate_sites.bed

    mosdepth --fast-mode --no-per-base -t ${task.cpus} \
        --by candidate_sites.bed sample ${bam}

    zcat sample.regions.bed.gz \
        | awk -v OFS='\t' -v min=${params.min_target_depth} '\$4 >= min {print \$1, \$2, \$3}' \
        > ${output_path}
    """
}

// Unions the per-sample kept-site BEDs from computeTargetDepth and filters
// the sorted target VCF down to that union.
//
// A site is kept if ANY sample clears the depth floor, not only if every
// sample does. Both the sharded and monolithic cellSNP paths below genotype
// every sample against one shared target list, so this must not drop a site
// that a later sample actually has coverage at just because an earlier
// sample lacked it -- per-sample-aware filtering still happens downstream via
// the per-SNP minMAF/minCOUNT filters in cellsnpDonorFilters. This step only
// removes sites that are dead weight for every sample in the run.
process filterTargetsByDepth {
    label 'small'
    publishDir "${params.output_dir}/cell_snp/",
        mode: 'copy',
        enabled: params.publish_cell_snp
    tag "Filter targets by depth"

    input:
    path depth_beds
    path snp_annotation

    output:
    path output_path

    script:
    output_path = "depth_filtered_snp_annotation.vcf.gz"
    """
    cat ${depth_beds} | sort -k1,1 -k2,2n -u > kept_sites.bed
    bcftools view -R kept_sites.bed ${snp_annotation} -Oz -o ${output_path}
    """
}

// Render a `path` input as cellsnp-lite's comma-separated -s argument.
//
// Do NOT call .join(',') on the input directly. A path input arrives as a bare
// java.nio.file.Path when it holds one file and only becomes a List when it
// holds several -- and Path is itself Iterable, over its *name elements*. So
// joining a Path silently yields "sub,dir,x.bam" for any staged path with
// directory components, and "" for an empty Path, which getopt then misreads as
// -s consuming the following flag. Normalise to a list of strings first.
//
// Shared by the monolithic and sharded processes so the two cannot drift apart.
def cellsnpBamArg(bam_paths) {
    def paths = (bam_paths instanceof Collection) ? bam_paths.toList() : [bam_paths]
    def names = paths.findAll { it != null }.collect { it.toString() }
    if (!names || names.any { it.isEmpty() }) {
        error "cellsnp-lite needs at least one BAM for -s, got: ${bam_paths}"
    }
    return names.join(",")
}

process runCellSNPGenotype {
    label 'large'
    publishDir "${params.output_dir}/cell_snp/",
        mode: 'copy',
        enabled: params.publish_cell_snp
    tag "CellSNP-genotype ${suffix}"

    input:
    tuple path(bam_paths), path(path_indices), path(barcode_path), val(suffix), val(n_donors), path(snp_annotation)

    output:
    tuple path(output_path), val(suffix)

    script:
    output_path = "cellsnp_" + suffix
    def bam_arg = cellsnpBamArg(bam_paths)
    def filters = cellsnpDonorFilters(n_donors)
    """
    cellsnp-lite -s "${bam_arg}" \
        -b ${barcode_path} \
        -O ${output_path} \
        -T ${snp_annotation} \
        -p ${task.cpus} \
        --minMAF ${filters.minMAF} \
        --minCOUNT ${filters.minCOUNT} \
        --gzip --genotype
    """
}

process splitCellSNPTargets {
    label 'small'
    tag "CellSNP-split-targets"

    input:
    tuple path(snp_annotation), val(n_chunks)

    output:
    path "part_*.vcf"

    script:
    """
    split_cellsnp_targets.py ${snp_annotation} ${n_chunks} part_
    """
}

process runCellSNPGenotypeChunk {
    label 'cellsnp_chunk'
    array params.cellsnp_chunks
    scratch true
    tag "CellSNP-genotype ${suffix} chunk ${chunk_index}"

    input:
    tuple path(bam_paths), path(path_indices), path(barcode_path), val(suffix), val(n_donors), val(chunk_index), path(snp_chunk)

    output:
    tuple val(suffix), val(chunk_index), path(output_path)

    script:
    output_path = "cellsnp_" + suffix + "_chunk_" + chunk_index
    def bam_arg = cellsnpBamArg(bam_paths)
    def filters = cellsnpDonorFilters(n_donors)
    """
    cellsnp-lite -s "${bam_arg}" \
        -b ${barcode_path} \
        -O ${output_path} \
        -T ${snp_chunk} \
        -p ${task.cpus} \
        --minMAF ${filters.minMAF} \
        --minCOUNT ${filters.minCOUNT} \
        --gzip --genotype
    """
}

process mergeCellSNP {
    label 'small'
    publishDir "${params.output_dir}/cell_snp/",
        mode: 'copy',
        enabled: params.publish_cell_snp
    tag "CellSNP-merge ${suffix}"

    input:
    tuple val(suffix), path(chunk_dirs), path(contig_order)

    output:
    tuple path(output_path), val(suffix)

    script:
    output_path = "cellsnp_" + suffix
    """
    merge_cellsnp.py ${output_path} ${contig_order} ${chunk_dirs}
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
