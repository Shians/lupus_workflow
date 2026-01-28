// Import modules
include { validateParams } from './modules/validation/main.nf'
include { createChannelFromSampleSheet; parseVireoSampleSheet } from './modules/sample_sheet/main.nf'
include { buildMinimapIndexGenome; buildMinimapIndexTranscriptome } from './modules/indexing/main.nf'
include { bamToFastq; splitFastqChunks } from './modules/preprocessing/main.nf'
include { flexiplexGetBarcodeCandidates; mergeFlexiplexBarcodes; flexiplexTagFastq } from './modules/barcode_detection/main.nf'
include { alignMinimap2Spliced; alignMinimap2TranscriptomeUnsorted } from './modules/alignment/main.nf'
include { catTranscriptAlignedBams; mergeSpliceAlignedBams; combineMergedSplicedBams; sortBamByName } from './modules/postprocessing/main.nf'
include { runCellSNPGenotype; runVireoDemultiplex } from './modules/demultiplexing/main.nf'
include { runOarfish } from './modules/quantification/main.nf'
include { runBcl2Fastq } from './modules/utilities/main.nf'

def helpMessage() {
    log.info """
    ===================================
    Lupus Workflow - Single Cell Analysis Pipeline
    ===================================

    Usage:
      nextflow run main.nf [options]

    Required Parameters:
      --sample_sheet PATH              Path to sample sheet TSV file
      --reference_genome PATH          Path to reference genome FASTA file
      --reference_transcriptome PATH   Path to reference transcriptome FASTA file
      --canonical_barcode_list PATH    Path to canonical barcode list file

    Optional Parameters:
      --output_dir PATH                Output directory (default: output)
      --snp_annotation PATH            Path to VCF file with SNP positions for CellSNP genotyping
      --vireo_sample_sheet PATH        Path to Vireo sample sheet TSV file specifying donor counts per sample

    Sample Sheet Format:
      The sample sheet must be a tab-separated (TSV) file with the following columns:
        - sample_id: Unique identifier for each sample (can include flowcell ID)
        - bam_dir: Directory containing BAM files for this sample

      The workflow will:
        1. Find all *.bam files in each directory
        2. Process each BAM file in parallel (BAM → FASTQ → alignment)
        3. Aggregate results by sample_id after alignment

      Example (columns separated by tabs):
        sample_id\tbam_dir
        Sample1_FC001\t/data/flowcell_FC001/sample1/
        Sample1_FC002\t/data/flowcell_FC002/sample1/
        Sample2_FC001\t/data/flowcell_FC001/sample2/
        Sample2_FC002\t/data/flowcell_FC002/sample2/

      This allows maximum parallelization: each BAM file is processed as a separate
      task, ideal for single-cell ONT data with hundreds of BAMs per flowcell.

    Genotyping and Demultiplexing (Optional):
      CellSNP genotyping runs when --snp_annotation is provided to generate variant call data.
      Vireo demultiplexing is optional and only runs when both --snp_annotation and
      --vireo_sample_sheet are provided.

    Vireo Sample Sheet Format (Optional):
      The Vireo sample sheet is a TSV file with the following columns:
        - sample_id: Must match sample_id values in the main sample sheet
        - n_donors: Number of donors expected in each sample (positive integer)

      Example (columns separated by tabs):
        sample_id\tn_donors
        Sample1\t3
        Sample2\t2

      The sample_id should match the base sample identifier without flowcell suffix.
      If a sample is not listed, it will default to 2 donors.

    Alignment Options:
      --alignment.bam_parts INT        Number of BAM parts for processing (default: 32)

    Barcode Detection Options:
      --barcode_detection.min_barcode_count INT   Minimum read count for a barcode to be included (default: 500)

    Publishing Options:
      Control which intermediate and final results are saved:

      Intermediate Files:
        --publish.fastq BOOL                   Save BAM to FASTQ converted files (default: false)
        --publish.barcoded_fastq BOOL          Save barcoded FASTQ files (default: true)
        --publish.flexiplex_candidates BOOL    Save individual barcode candidates (default: false)
        --publish.flexiplex_merged BOOL        Save merged barcode list (default: true)
        --publish.flexiplex_logs BOOL          Save Flexiplex log files (default: true)
        --publish.splice_aligned BOOL          Save splice-aligned BAM files (default: true)

      Final Results:
        --publish.oarfish BOOL                 Save Oarfish quantification results (default: true)
        --publish.cell_snp BOOL                Save CellSNP genotyping results (default: true)
        --publish.vireo BOOL                   Save Vireo demultiplexing results (default: true)

      Utility Outputs:
        --publish.bcl2fastq BOOL               Save BCL to FASTQ conversion results (default: false)

    Execution Profiles:
      -profile conda                   Use Conda for dependency management
      -profile docker                  Use Docker for containerization
      -profile singularity             Use Singularity for containerization

    Other Options:
      --help                           Display this help message

    Example:
      nextflow run main.nf \\
        --sample_sheet samples.tsv \\
        --reference_genome /path/to/genome.fa \\
        --reference_transcriptome /path/to/transcriptome.fa \\
        --canonical_barcode_list /path/to/barcodes.txt \\
        --output_dir results \\
        -profile conda
    """
}

workflow {
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Validate input parameters before starting the workflow
    validateParams()
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
        .map { sample, fastq -> tuple(sample, fastq, params.alignment.bam_parts) }
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
