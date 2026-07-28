include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

// Import modules
include { createChannelFromSampleSheet; parseVireoSampleSheet } from './modules/sample_sheet/main.nf'
include { buildMinimapIndexGenome; buildMinimapIndexTranscriptome } from './modules/indexing/main.nf'
include { bamToFastq; splitFastqChunks } from './modules/preprocessing/main.nf'
include { flexiplexGetBarcodeCandidates; mergeFlexiplexBarcodes; flexiplexTagFastq } from './modules/barcode_detection/main.nf'
include { alignMinimap2Spliced; alignMinimap2TranscriptomeUnsorted } from './modules/alignment/main.nf'
include { catTranscriptAlignedBams; mergeSpliceAlignedBams; sortBamByCellBarcode } from './modules/postprocessing/main.nf'
include { indexBam } from './modules/deduplication/main.nf'
include { sortIndexBam } from './modules/deduplication/main.nf'
include { nailpolishDedup } from './modules/deduplication/main.nf'
include { umitoolsDedup as umitoolsDedupGenome } from './modules/deduplication/main.nf'
include { umitoolsDedup as umitoolsDedupTranscriptome } from './modules/deduplication/main.nf'
include { sortSNPAnnotation; computeTargetDepth; filterTargetsByDepth; runCellSNPGenotype; splitCellSNPTargets; runCellSNPGenotypeChunk; mergeCellSNP; runVireoDemultiplex } from './modules/demultiplexing/main.nf'
include { runOarfish } from './modules/quantification/main.nf'
include { craminoStats; samtoolsFlagstat; mosdepthCoverage; readCountSummary; multiQC } from './modules/qc/main.nf'
include { countReads as countRawReads } from './modules/qc/main.nf'
include { countReads as countBarcodeTaggedReads } from './modules/qc/main.nf'
include { countReads as countGenomeAlignedReads } from './modules/qc/main.nf'
include { countReads as countGenomeDedupReads } from './modules/qc/main.nf'

workflow {
    // Validate parameters against schema and print summary
    validateParameters()
    log.info paramsSummaryLog(workflow)

    ref_genome_path = channel.fromPath(params.reference_genome)
    transcriptome_path = channel.fromPath(params.reference_transcriptome)
    canonical_bc_list = channel.fromPath(params.canonical_barcode_list)

    ref_genome_index = buildMinimapIndexGenome(ref_genome_path)
    transcriptome_index = buildMinimapIndexTranscriptome(transcriptome_path)

    // Create FASTQ files from input BAM files
    // Each BAM file is processed individually for parallelization
    // Parse TSV and expand directories to individual BAM files
    // Channel contains: [sample_id, bam_file] for each BAM file
    bam_channel = createChannelFromSampleSheet(params.sample_sheet)

    // QC: Cramino stats and raw read counts run in parallel with BAM conversion
    cramino_out = craminoStats(bam_channel)
    raw_counts_ch = countRawReads(bam_channel.map { sample_id, bam -> tuple(sample_id, bam, 'input_bam', 'PRIMARY') })

    // Convert each BAM to FASTQ in parallel
    untagged_fastq_files = bamToFastq(bam_channel)

    // Barcode detection - process each FASTQ individually
    flexiplex_candidate_bc = flexiplexGetBarcodeCandidates(untagged_fastq_files)

    // Merge candidate barcodes with canonical list for each sample
    // Group all barcode count files by sample_id before merging
    merged_bc_result = flexiplex_candidate_bc
        .groupTuple()
        .combine(canonical_bc_list)
        | mergeFlexiplexBarcodes
    flexiplex_bc = merged_bc_result.barcodes

    // Tag FASTQ files with detected barcodes
    tagged_fastq_files = untagged_fastq_files
        .combine(flexiplex_bc, by: 0)
        | flexiplexTagFastq

    // Deduplication strategy (see params.dedup_method):
    //   'nailpolish' - molecular (CB+UMI) consensus dedup at the FASTQ level, before
    //                  alignment. Each molecule aligns once, preserving the full
    //                  multimapping set for oarfish, and serves both the genome
    //                  (CellSNP) and transcriptome (oarfish) paths from one step.
    //   'umitools'   - legacy post-alignment position-based dedup, run separately on
    //                  the genome and transcriptome BAMs (kept for validation).
    if (params.dedup_method == 'nailpolish') {
        reads_for_alignment = nailpolishDedup(tagged_fastq_files.fastq).fastq
    } else {
        reads_for_alignment = tagged_fastq_files.fastq
    }

    // Split FASTQ into chunks for parallel alignment
    fastq_chunks = reads_for_alignment
        .map { sample, fastq -> tuple(sample, fastq, params.bam_parts) }
        | splitFastqChunks
        | transpose

    // Genome alignment
    merged_spliced_bams = fastq_chunks.combine(ref_genome_index)
        | alignMinimap2Spliced
        | groupTuple
        | mergeSpliceAlignedBams

    // QC: genome alignment stats and read count
    flagstat_out = samtoolsFlagstat(merged_spliced_bams)
    mosdepth_out = mosdepthCoverage(merged_spliced_bams)
    // In nailpolish mode the aligned BAM is already deduplicated (dedup happened on the
    // FASTQ), so this first count reflects deduplicated molecules; in umitools mode it
    // reflects reads that entered alignment, before the post-alignment dedup below.
    aligned_stage = params.dedup_method == 'nailpolish' ? 'after_dedup_aligned' : 'after_barcode_tagging'
    tagged_counts_ch = countBarcodeTaggedReads(merged_spliced_bams.map { sample_id, bam, _bai -> tuple(sample_id, bam, aligned_stage, 'PRIMARY') })
    genome_counts_ch = countGenomeAlignedReads(merged_spliced_bams.map { sample_id, bam, _bai -> tuple(sample_id, bam, 'after_genome_alignment', 'PRIMARY_MAPPED') })

    // Genome dedup + indexed BAM for CellSNP.
    if (params.dedup_method == 'nailpolish') {
        // Reads were already deduplicated at the FASTQ level; the merged spliced BAM is
        // already coordinate-sorted and indexed by mergeSpliceAlignedBams.
        genome_bam_indexed = merged_spliced_bams
        dedup_genome_counts_ch = channel.empty()
    } else {
        dedup_bam_genome = umitoolsDedupGenome(merged_spliced_bams)
        dedup_genome_counts_ch = countGenomeDedupReads(dedup_bam_genome.bam.map { sample_id, bam -> tuple(sample_id, bam, 'after_genome_dedup', 'PRIMARY_MAPPED') })
        genome_bam_indexed = indexBam(dedup_bam_genome.bam).bam_bai
    }

    // Transcript alignment
    transcript_aligned = fastq_chunks.combine(transcriptome_index)
        | alignMinimap2TranscriptomeUnsorted
        | groupTuple
        | catTranscriptAlignedBams

    // Transcript dedup: nailpolish already deduplicated pre-alignment, so feed the
    // concatenated BAM straight through; umitools dedups the sorted transcriptome BAM.
    if (params.dedup_method == 'nailpolish') {
        transcript_bam = transcript_aligned
    } else {
        sorted_indexed_transcript = sortIndexBam(transcript_aligned)
        transcript_bam = umitoolsDedupTranscriptome(sorted_indexed_transcript).bam
    }

    // oarfish requires the BAM collated by cell barcode
    transcript_bam
        | sortBamByCellBarcode
        | runOarfish

    // Read tracking summary
    all_counts_ch = raw_counts_ch
        .mix(tagged_counts_ch)
        .mix(genome_counts_ch)
        .mix(dedup_genome_counts_ch)
        .groupTuple()
    readCountSummary(all_counts_ch)

    // MultiQC: collect cramino, flagstat, and mosdepth outputs
    qc_inputs = cramino_out.map { _id, json -> json }
        .mix(flagstat_out.map { _id, txt -> txt })
        .mix(mosdepth_out.map { _id, summary, dist -> [summary, dist] })
        .flatten()
        .collect()
    multiQC(qc_inputs)

    // Genotyping and Demultiplexing
    if (params.snp_annotation) {
        // Sorted once here, ahead of the monolithic/sharded branch below, so
        // both paths genotype the same coordinate-ordered site list. The
        // contig order is carried alongside as a value channel: mergeCellSNP
        // needs it to rebuild that same ordering, and every merge task reads
        // the same copy.
        sorted_snp_annotation = sortSNPAnnotation(channel.fromPath(params.snp_annotation))
        sorted_snp_annotation_vcf = sorted_snp_annotation.vcf
        contig_order_path = sorted_snp_annotation.contig_order.first()

        // Donor counts drive both the Vireo -N argument and the cellsnp-lite
        // site filters (see cellsnpDonorFilters in modules/demultiplexing), so
        // resolve them once here instead of only inside the Vireo branch below.
        // Samples missing from the sheet, or a run with no sheet at all, fall
        // back to two donors -- the previous Vireo default.
        default_n_donors = 2
        vireo_donors = params.vireo_sample_sheet
            ? parseVireoSampleSheet(params.vireo_sample_sheet)
            : [:]

        // Depth-filter the sorted target list once here, ahead of both the
        // monolithic and sharded branches below, so every genotyping job --
        // whichever branch runs -- genotypes the exact same filtered VCF.
        // See computeTargetDepth/filterTargetsByDepth in
        // modules/demultiplexing for why this is cheap relative to the
        // cellSNP genotyping it is trimming.
        if (params.filter_targets_by_depth) {
            depth_beds = genome_bam_indexed
                .combine(sorted_snp_annotation_vcf)
                .map { sample_id, bam, bai, vcf -> tuple(bam, bai, sample_id, vcf) }
                | computeTargetDepth
                | collect

            snp_annotation_path = filterTargetsByDepth(depth_beds, sorted_snp_annotation_vcf)
        } else {
            snp_annotation_path = sorted_snp_annotation_vcf
        }

        // Prepare per-sample (bam, bai, barcode, sample_id, n_donors) tuples
        // for CellSNP. genome_bam_indexed is (sample_id, bam, bai) in both
        // dedup modes.
        cellsnp_bam_bc = genome_bam_indexed
            .join(flexiplex_bc)
            .map { sample_id, bam, bai, barcode_file ->
                tuple(bam, bai, barcode_file, sample_id, vireo_donors[sample_id] ?: default_n_donors)
            }

        // CellSNP genotyping: by default (params.cellsnp_sharded) split the
        // target VCF into balanced chunks, genotype each as an independent
        // right-sized array job, and merge back into a single cellSNP
        // directory (identical to a single run, so Vireo consumes it
        // unchanged). Set cellsnp_sharded = false for one monolithic job.
        // See modules/demultiplexing.
        if (params.cellsnp_sharded) {
            // One tuple per chunk VCF, tagged with its zero-padded index taken
            // from the "part_<idx>.vcf" filename so ordering is deterministic.
            target_chunks = splitCellSNPTargets(
                    snp_annotation_path.map { vcf -> tuple(vcf, params.cellsnp_chunks) }
                )
                .flatMap { chunk_files -> (chunk_files instanceof List) ? chunk_files : [chunk_files] }
                .map { chunk_vcf ->
                    def chunk_index = (chunk_vcf.name =~ /part_(\d+)\.vcf/)[0][1]
                    tuple(chunk_index, chunk_vcf)
                }

            // Cartesian product: every sample against every target chunk.
            cellsnp_chunk_input = cellsnp_bam_bc
                .combine(target_chunks)
                .map { bam, bai, barcode_file, sample_id, n_donors, chunk_index, chunk_vcf ->
                    tuple(bam, bai, barcode_file, sample_id, n_donors, chunk_index, chunk_vcf)
                }

            cellsnp_results = runCellSNPGenotypeChunk(cellsnp_chunk_input)
                .map { suffix, _chunk_index, chunk_dir -> tuple(suffix, chunk_dir) }
                .groupTuple()
                .combine(contig_order_path)
                | mergeCellSNP
        } else {
            cellsnp_input = cellsnp_bam_bc.combine(snp_annotation_path)
            cellsnp_results = runCellSNPGenotype(cellsnp_input)
        }

        // Conditionally run Vireo demultiplexing
        if (params.vireo_sample_sheet) {
            cellsnp_results
                .map { cellsnp_path, sample_id ->
                    tuple(cellsnp_path, sample_id, vireo_donors[sample_id] ?: default_n_donors)
                }
                | runVireoDemultiplex
        }
    }
}
