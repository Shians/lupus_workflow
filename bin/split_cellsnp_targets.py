#!/usr/bin/env python3
"""Split a cellSNP target SNP VCF into N chunks with ~equal record counts.

Each chunk carries the full VCF header (all ``##`` meta lines plus the
``#CHROM`` line) so that cellsnp-lite accepts it as a ``-T`` target file.

Work in cellsnp-lite Mode 1a scales with the number of SNPs times depth, so
chunks are balanced by *record count*, not by base-pair span. Records are
distributed round-robin across chunks: this interleaves chromosomes so that a
single heavy chromosome (e.g. chr1) does not concentrate in one chunk, giving
better wall-clock balance than contiguous slicing while keeping counts equal.

Usage:
    split_cellsnp_targets.py <input.vcf[.gz]> <n_chunks> <out_prefix>

Emits files ``<out_prefix>000.vcf``, ``<out_prefix>001.vcf``, ... one per
non-empty chunk. Output is plain (uncompressed) VCF, which cellsnp-lite reads
directly with no index required.
"""

import gzip
import sys


def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)

    input_vcf = sys.argv[1]
    requested_chunks = int(sys.argv[2])
    out_prefix = sys.argv[3]

    if requested_chunks < 1:
        sys.exit(f"n_chunks must be >= 1, got {requested_chunks}")

    header_lines = []
    record_count = 0

    # First pass: read the header and count records so we can cap the number of
    # chunks (never emit an empty target VCF, which cellsnp-lite may reject).
    with open_maybe_gz(input_vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                record_count += 1

    if record_count == 0:
        sys.exit(f"No variant records found in {input_vcf}")

    if not any(line.startswith("#CHROM") for line in header_lines):
        sys.exit(f"No #CHROM line found in header of {input_vcf}")

    n_chunks = min(requested_chunks, record_count)

    # Open all chunk handles and write the shared header to each.
    handles = []
    for chunk_index in range(n_chunks):
        path = f"{out_prefix}{chunk_index:03d}.vcf"
        handle = open(path, "wt")
        handle.writelines(header_lines)
        handles.append(handle)

    # Second pass: distribute records round-robin across the chunk handles.
    record_index = 0
    with open_maybe_gz(input_vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            handles[record_index % n_chunks].write(line)
            record_index += 1

    for handle in handles:
        handle.close()

    print(
        f"Split {record_count} SNP records into {n_chunks} chunks "
        f"(requested {requested_chunks})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
