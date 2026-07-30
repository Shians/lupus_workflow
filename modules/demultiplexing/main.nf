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
// sets, and it is the precondition for merging sharded chunks back into
// coordinate order: splitCellSNPTargets cuts this sorted list into contiguous
// ranges, so the chunks come back disjoint and already ordered and mergeCellSNP
// can concatenate them rather than reorder them.
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
//
// NOT --fast-mode. mosdepth's help for that flag is "dont look at internal
// cigar operations", which makes it count a read across its whole POS->end
// span, introns included. On spliced RNA-seq that is a wildly different number
// from base coverage: a single long read bridging a 100 kb intron would credit
// every site under that intron with depth 1. Filtering on span depth means
// sites whose only "coverage" is introns arching over them survive
// min_target_depth, inflating the target list with sites that have no usable
// reads -- and those sites then cost a full per-cell pileup before
// cellsnp-lite's own minCOUNT discards them, which is the exact waste this
// sieve exists to prevent. The floor is about usable bases, so it needs the
// CIGAR-aware count.
//
// Dropping --fast-mode costs runtime in this process, which is cheap relative
// to the per-cell cellSNP pass it is protecting.
//
// This depth feeds the site filter and nothing else. A previous revision also
// used it to weight the shard split; that was removed once a full array showed
// shard runtime is uncorrelated with read depth (+0.04) and tracks site count
// instead -- see split_cellsnp_targets.py, which now balances on counts and
// needs no depth input at all.
process computeTargetDepth {
    label 'medium'
    tag "Target-depth ${sample_id}"

    input:
    tuple path(bam), path(bai), val(sample_id), path(snp_annotation)

    output:
    path output_path

    script:
    output_path = "depth_${sample_id}.bed.gz"
    """
    set -o pipefail

    bcftools query -f '%CHROM\\t%POS\\n' ${snp_annotation} \
        | awk -v OFS='\t' '{print \$1, \$2-1, \$2}' > candidate_sites.bed

    mosdepth --no-per-base -t ${task.cpus} \
        --by candidate_sites.bed sample ${bam}

    # mosdepth groups its output by the BAM header's contig order, which is not
    # required to match the target VCF's. aggregate_target_depth.py reads the
    # per-sample tables row-for-row, so they must agree on row order; it fails
    # loudly on a mismatch rather than mispairing sites, but it fails late, once
    # every sample's depth pass has already run. Check the precondition here,
    # where it is one cheap pass and the error names the real cause.
    paste candidate_sites.bed <(zcat sample.regions.bed.gz) \
        | awk '\$1 != \$4 || \$2 != \$5 {
                print "mosdepth emitted regions in a different order than " \
                      "candidate_sites.bed at line " NR ", so the BAM header " \
                      "contig order differs from the target VCF" > "/dev/stderr"
                exit 1
            }'

    # Passed on whole rather than filtered here: filterTargetsByDepth needs
    # every sample's depth at every candidate site, both to apply the
    # keep-if-any rule and to take the per-site maximum that weights the shards.
    mv sample.regions.bed.gz ${output_path}
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
    set -o pipefail

    aggregate_target_depth.py ${params.min_target_depth} ${params.min_target_sites} \
        kept_sites.bed ${depth_beds}

    # -T, not -R: -R index-jumps and so requires a .csi/.tbi that
    # sortSNPAnnotation does not produce, while -T streams the whole VCF and
    # needs no index. Streaming is the right access pattern anyway -- the kept
    # set is a large fraction of a genome-wide list, not a handful of regions --
    # and it emits records in VCF order rather than in kept_sites.bed order,
    # which is what keeps the output coordinate-sorted for splitCellSNPTargets
    # and mergeCellSNP. The .bed suffix is load-bearing: it is how bcftools
    # knows to read the file as 0-based half-open rather than 1-based.
    bcftools view -T kept_sites.bed ${snp_annotation} -Oz -o ${output_path}
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

// NOT --genotype. That flag is what makes cellsnp-lite emit
// cellSNP.cells.vcf.gz, a per-cell GT/AD/DP/PL matrix laid out as one VCF column
// per barcode. At 41k cells it is pathologically compressible -- 33 MB on disk
// per shard, 5.5 GB inflated (167x), because nearly every cell carries a
// `.:.:.:.` placeholder on nearly every row -- and no consumer opens it:
// vireo's read_cellSNP and snplet::import_cellsnp both read only base.vcf.gz,
// the AD/DP/OTH matrices and samples.tsv. The `-t GT` on runVireoDemultiplex
// names the tag of a donor genotype *file* (-d), which is not passed, so it
// does not make cells.vcf load-bearing either.
//
// The AD/DP/OTH counting matrices, which are the actual input to vireo, are
// emitted regardless of this flag; --genotype only adds the per-cell VCF and
// the genotype-likelihood work behind it. Dropping it cost ~64 min of merging
// and shrinks every shard's own output.
//
// Restore it only to produce cells.vcf as a donor-genotype reference for a
// future `vireo -d` run -- merge_cellsnp.py handles its presence or absence.
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
        --maxDEPTH ${params.cellsnp_max_depth} \
        --maxPILEUP ${params.cellsnp_max_pileup} \
        --gzip
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
        --maxDEPTH ${params.cellsnp_max_depth} \
        --maxPILEUP ${params.cellsnp_max_pileup} \
        --gzip
    """
}

// Reassemble the sharded genotyping into one cellSNP directory for vireo.
//
// The chunks are contiguous ranges of the sorted target VCF, so this is a
// concatenation with a row offset per chunk, not a k-way merge -- see
// merge_cellsnp.py, which verifies that contiguity rather than assuming it.
//
// The cores are for compression, not for the merge, which is inherently serial.
// merge_cellsnp.py hands --threads to bgzip (preferred, since it restores the
// BGZF blocking that cellsnp-lite wrote and plain gzip would discard) or pigz,
// falling back to single-core zlib when neither is on PATH -- correct either
// way, only the runtime changes.
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
    merge_cellsnp.py ${output_path} ${contig_order} ${chunk_dirs} \
        --threads ${task.cpus}
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
