// Import modules
include { validateParams } from './modules/validation/main.nf'
include { buildMinimapIndexGenome; buildMinimapIndexTranscriptome } from './modules/indexing/main.nf'
include { bamToFastq; splitFastqChunks; retagFastq } from './modules/preprocessing/main.nf'
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
      --input_bam PATH                 Path to input BAM file
      --reference_genome PATH          Path to reference genome FASTA file
      --reference_transcriptome PATH   Path to reference transcriptome FASTA file
      --canonical_barcode_list PATH    Path to canonical barcode list file

    Optional Parameters:
      --output_dir PATH                Output directory (default: output)

    Alignment Options:
      --alignment.bam_parts INT        Number of BAM parts for processing (default: 32)

    Publishing Options:
      Control which intermediate and final results are saved:

      Intermediate Files:
        --publish.fastq BOOL                   Save BAM to FASTQ converted files (default: false)
        --publish.retagged_fastq BOOL          Save retagged FASTQ files (default: true)
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
    ref_genome_path = Channel.fromPath(params.reference_genome)
    transcriptome_path = Channel.fromPath(params.reference_transcriptome)
    canonical_bc_list = Channel.fromPath(params.canonical_barcode_list)

    ref_genome_index = buildMinimapIndexGenome(ref_genome_path)
    transcriptome_index = buildMinimapIndexTranscriptome(transcriptome_path)

    untagged_fastq_files = channel.from("Untreated1")
        .combine(channel.fromPath("data/fastq/*.fastq.gz"))

    flexiplex_candidate_bc = flexiplexGetBarcodeCandidates(untagged_fastq_files)
        .map{ sample, bc_file -> bc_file }
        .collect()

    flexiplex_bc = mergeFlexiplexBarcodes(flexiplex_candidate_bc, canonical_bc_list)

    tagged_fastq_files = untagged_fastq_files
        .combine(flexiplex_bc)
        | flexiplexTagFastq

    fastq_chunks = tagged_fastq_files.fastq
        .map { sample, fastq -> tuple(sample, fastq, 10) }
        | splitFastqChunks
        | transpose

    retagged_fastq = retagFastq(fastq_chunks)

    merged_spliced_bams = retagged_fastq.combine(ref_genome_index)
        | alignMinimap2Spliced
        | groupTuple
        | mergeSpliceAlignedBams

    // Transcript quantification
    // retagged_fastq.combine(transcriptome_index)
    //     | alignMinimap2TranscriptomeUnsorted
    //     | groupTuple
    //     | catTranscriptAlignedBams
    //     | sortBamByName
    //     | runOarfish

    // runCellSNPGenotype(cellsnp_input) |
    //     runVireoDemultiplex
}
