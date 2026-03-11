include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

// Import modules
include { createChannelFromSampleSheet; parseVireoSampleSheet } from './modules/sample_sheet/main.nf'
include { buildMinimapIndexGenome; buildMinimapIndexTranscriptome } from './modules/indexing/main.nf'
include { bamToFastq; splitFastqChunks } from './modules/preprocessing/main.nf'
include { flexiplexGetBarcodeCandidates; mergeFlexiplexBarcodes; flexiplexTagFastq } from './modules/barcode_detection/main.nf'
include { alignMinimap2Spliced; alignMinimap2TranscriptomeUnsorted } from './modules/alignment/main.nf'
include { catTranscriptAlignedBams; mergeSpliceAlignedBams; combineMergedSplicedBams; sortBamByName } from './modules/postprocessing/main.nf'
include { runCellSNPGenotype; runVireoDemultiplex } from './modules/demultiplexing/main.nf'
include { runOarfish } from './modules/quantification/main.nf'

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

    // Convert each BAM to FASTQ in parallel
    untagged_fastq_files = bamToFastq(bam_channel)

    // Barcode detection - process each FASTQ individually
    flexiplex_candidate_bc = flexiplexGetBarcodeCandidates(untagged_fastq_files)

    // Merge candidate barcodes with canonical list for each sample
    // Group all barcode count files by sample_id before merging
    flexiplex_bc = flexiplex_candidate_bc
        .groupTuple()
        .combine(canonical_bc_list)
        | mergeFlexiplexBarcodes

    // Tag FASTQ files with detected barcodes
    tagged_fastq_files = untagged_fastq_files
        .join(flexiplex_bc)
        | flexiplexTagFastq

    // Split FASTQ into chunks for parallel alignment
    fastq_chunks = tagged_fastq_files.fastq
        .map { sample, fastq -> tuple(sample, fastq, params.bam_parts) }
        | splitFastqChunks
        | transpose

    // Genome alignment
    merged_spliced_bams = fastq_chunks.combine(ref_genome_index)
        | alignMinimap2Spliced
        | groupTuple
        | mergeSpliceAlignedBams

    // Transcript quantification
    fastq_chunks.combine(transcriptome_index)
        | alignMinimap2TranscriptomeUnsorted
        | groupTuple
        | catTranscriptAlignedBams
        | sortBamByName
        | runOarfish

    // Genotyping and Demultiplexing
    if (params.snp_annotation) {
        snp_annotation_path = channel.fromPath(params.snp_annotation)

        // Prepare input for CellSNP: combine BAMs with their barcode lists
        cellsnp_input = merged_spliced_bams
            .join(flexiplex_bc)
            .map { sample_id, bam, bai, barcode_file ->
                tuple([bam], [bai], barcode_file, sample_id)
            }
            .combine(snp_annotation_path)

        // Run CellSNP genotyping (always runs when snp_annotation is provided)
        cellsnp_results = runCellSNPGenotype(cellsnp_input)

        // Conditionally run Vireo demultiplexing
        if (params.vireo_sample_sheet) {
            vireo_donors = parseVireoSampleSheet(params.vireo_sample_sheet)

            cellsnp_results
                .map { cellsnp_path, sample_id ->
                    tuple(cellsnp_path, sample_id, vireo_donors[sample_id] ?: 2)
                }
                | runVireoDemultiplex
        }
    }
}
