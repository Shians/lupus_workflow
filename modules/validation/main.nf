// Input validation functions for pipeline parameters

// Helper function to validate a required file parameter
def validateFileParam(errors, paramValue, paramName, displayName) {
    if (!paramValue) {
        errors << "ERROR: ${paramName} is not defined"
    } else if (paramValue.trim().isEmpty()) {
        errors << "ERROR: ${paramName} cannot be empty"
    } else {
        def fileObj = file(paramValue)
        if (!fileObj.exists()) {
            errors << "ERROR: ${displayName} does not exist: ${paramValue}"
        } else if (fileObj.isEmpty()) {
            errors << "ERROR: ${displayName} is empty: ${paramValue}"
        } else if (!fileObj.isFile()) {
            errors << "ERROR: ${displayName} path is not a file: ${paramValue}"
        }
    }
}

// Helper function to validate an optional file parameter
def validateOptionalFileParam(errors, paramValue, paramName, displayName) {
    if (paramValue && !paramValue.trim().isEmpty()) {
        def fileObj = file(paramValue)
        if (!fileObj.exists()) {
            errors << "ERROR: ${displayName} does not exist: ${paramValue}"
        } else if (fileObj.isEmpty()) {
            errors << "ERROR: ${displayName} is empty: ${paramValue}"
        } else if (!fileObj.isFile()) {
            errors << "ERROR: ${displayName} path is not a file: ${paramValue}"
        }
    }
}

def validateParams() {
    def errors = []

    // Validate required file inputs
    validateFileParam(errors, params.sample_sheet, "params.sample_sheet", "Sample sheet file")
    validateFileParam(errors, params.reference_genome, "params.reference_genome", "Reference genome file")
    validateFileParam(errors, params.reference_transcriptome, "params.reference_transcriptome", "Reference transcriptome file")
    validateFileParam(errors, params.canonical_barcode_list, "params.canonical_barcode_list", "Canonical barcode list file")

    // Validate optional file inputs
    validateOptionalFileParam(errors, params.vireo_sample_sheet, "params.vireo_sample_sheet", "Vireo sample sheet file")
    validateOptionalFileParam(errors, params.snp_annotation, "params.snp_annotation", "SNP annotation VCF file")

    // Validate alignment parameters
    if (!params.alignment) {
        errors << "ERROR: params.alignment configuration is not defined"
    } else {
        // Validate bam_parts
        if (!params.alignment.bam_parts) {
            errors << "ERROR: params.alignment.bam_parts is not defined"
        } else if (!(params.alignment.bam_parts instanceof Integer)) {
            errors << "ERROR: params.alignment.bam_parts must be an integer, got: ${params.alignment.bam_parts}"
        } else if (params.alignment.bam_parts <= 0) {
            errors << "ERROR: params.alignment.bam_parts must be positive, got: ${params.alignment.bam_parts}"
        } else if (params.alignment.bam_parts > 1000) {
            errors << "WARNING: params.alignment.bam_parts is very high (${params.alignment.bam_parts}), this may cause performance issues"
        }
    }

    // Validate output_dir
    if (!params.output_dir) {
        errors << "ERROR: params.output_dir is not defined"
    } else if (params.output_dir.trim().isEmpty()) {
        errors << "ERROR: params.output_dir cannot be empty"
    }

    // Print all errors
    if (errors.size() > 0) {
        log.error "Parameter validation failed with ${errors.size()} error(s):"
        errors.each { error ->
            log.error "  ${error}"
        }
        System.exit(1)
    }

    // Print validation success message
    log.info "✓ Parameter validation successful"
    log.info "  - Output directory: ${params.output_dir}"
    log.info "  - Sample sheet: ${params.sample_sheet}"
    log.info "  - Reference genome: ${params.reference_genome}"
    log.info "  - Reference transcriptome: ${params.reference_transcriptome}"
    log.info "  - Canonical barcode list: ${params.canonical_barcode_list}"
    log.info "  - BAM parts: ${params.alignment.bam_parts}"
}
