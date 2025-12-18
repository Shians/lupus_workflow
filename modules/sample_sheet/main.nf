// Sample sheet parsing and validation functions

// Parse sample sheet TSV file
// Expected columns: sample_id, bam_dir
// This will expand each bam_dir to find all .bam files
def parseSampleSheet(sampleSheetPath) {
    def sampleSheet = []
    def lineNumber = 0
    def errors = []

    file(sampleSheetPath).withReader { reader ->
        def header = null
        reader.eachLine { line ->
            lineNumber += 1

            // Skip empty lines
            if (line.trim().isEmpty()) {
                return
            }

            // Parse header
            if (header == null) {
                header = line.split('\t').collect { it.trim() }

                // Validate required columns
                def requiredCols = ['sample_id', 'bam_dir']
                def missingCols = requiredCols.findAll { !header.contains(it) }
                if (missingCols) {
                    errors << "ERROR: Sample sheet missing required columns: ${missingCols.join(', ')}"
                    errors << "       Found columns: ${header.join(', ')}"
                    return
                }
                return
            }

            // Parse data rows
            def values = line.split('\t').collect { it.trim() }
            if (values.size() != header.size()) {
                errors << "ERROR: Line ${lineNumber} has ${values.size()} columns but header has ${header.size()} columns"
                return
            }

            def row = [header, values].transpose().collectEntries()

            // Validate sample_id
            if (!row.sample_id || row.sample_id.isEmpty()) {
                errors << "ERROR: Line ${lineNumber}: sample_id cannot be empty"
            }

            // Validate bam_dir exists
            if (!row.bam_dir || row.bam_dir.isEmpty()) {
                errors << "ERROR: Line ${lineNumber}: bam_dir cannot be empty"
            } else {
                def bamDir = file(row.bam_dir)
                if (!bamDir.exists()) {
                    errors << "ERROR: Line ${lineNumber}: BAM directory does not exist: ${row.bam_dir}"
                } else if (!bamDir.isDirectory()) {
                    errors << "ERROR: Line ${lineNumber}: BAM path is not a directory: ${row.bam_dir}"
                } else {
                    // Find all .bam files in the directory
                    def bamFiles = bamDir.listFiles().findAll { it.name.endsWith('.bam') }
                    if (bamFiles.isEmpty()) {
                        errors << "ERROR: Line ${lineNumber}: No BAM files found in directory: ${row.bam_dir}"
                    }
                }
            }

            sampleSheet << row
        }
    }

    if (errors.size() > 0) {
        log.error "Sample sheet validation failed with ${errors.size()} error(s):"
        errors.each { error ->
            log.error "  ${error}"
        }
        System.exit(1)
    }

    if (sampleSheet.size() == 0) {
        log.error "ERROR: Sample sheet contains no data rows"
        System.exit(1)
    }

    return sampleSheet
}

// Create a channel from sample sheet
// Returns a channel with tuples: [sample_id, bam_file_path]
// Each BAM file in each directory becomes a separate channel item for parallel processing
def createChannelFromSampleSheet(sampleSheetPath) {
    def sampleSheet = parseSampleSheet(sampleSheetPath)

    // Expand directories to individual BAM files
    def bamFilesList = []
    sampleSheet.each { row ->
        def bamDir = file(row.bam_dir)
        def bamFiles = bamDir.listFiles().findAll { it.name.endsWith('.bam') }
        bamFiles.each { bamFile ->
            bamFilesList << [sample_id: row.sample_id, bam_file: bamFile]
        }
    }

    // Log summary statistics
    def uniqueSamples = bamFilesList.collect { it.sample_id }.unique()
    def bamsBySample = bamFilesList.groupBy { it.sample_id }.collectEntries { k, v -> [k, v.size()] }

    log.info "Sample sheet loaded successfully:"
    log.info "  - Total BAM files: ${bamFilesList.size()}"
    log.info "  - Unique samples: ${uniqueSamples.size()}"
    log.info "  - Samples: ${uniqueSamples.join(', ')}"
    bamsBySample.each { sample, count ->
        log.info "    ${sample}: ${count} BAM files"
    }

    return Channel.fromList(
        bamFilesList.collect { row ->
            tuple(row.sample_id, file(row.bam_file))
        }
    )
}
