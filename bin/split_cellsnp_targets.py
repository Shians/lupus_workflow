#!/usr/bin/env python3
"""Split a cellSNP target SNP VCF into N contiguous chunks of equal record count.

Each chunk carries the full VCF header (all ``##`` meta lines plus the
``#CHROM`` line) so that cellsnp-lite accepts it as a ``-T`` target file.

Two properties matter, and they are independent:

**Data locality.** Chunks are *contiguous* runs of the coordinate-sorted target
list, so each shard reads one genomic interval of the BAM. An earlier
round-robin split gave every shard sites on every chromosome, so every shard
traversed the whole BAM: measured at ~17.7 TB of reads across 64 shards of a
276 GB BAM, against 1.0 TB once the split became contiguous.

**Balance.** Equal record counts, because runtime tracks the number of sites and
essentially nothing else. Measured over a full 64-shard array:

* ``corr(realtime, rchar) = +0.04`` -- runtime is unrelated to read volume, and
  equally unrelated to peak RSS (+0.03) and thread count (-0.13).
* The residual points at site count directly: one chunk ran 188 min on the
  *lowest* read volume in the array (5.7 GB) while another ran 5 min on 12 GB.
  Few-but-deep sites are fast; many-but-shallow sites are slow.

A previous revision weighted chunks by summed read depth instead. That is worse
than useless here: balancing sum-of-depth hands the many-shallow-sites chunks
proportionally more sites, which is the quantity that actually costs time. It
achieved 41% balance efficiency -- a 236 min critical path against a 98 min
perfectly-balanced ideal on 104 core-h of work. Do not reintroduce depth
weighting without first re-measuring the correlation above; if a cost model is
wanted, the term to fit is per-site, not per-read.

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


def scan_vcf(input_vcf):
    """Return ``(header_lines, record_count)`` for the target VCF."""
    header_lines = []
    record_count = 0
    with open_maybe_gz(input_vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                record_count += 1
    return header_lines, record_count


def chunk_sizes(record_count, n_chunks):
    """Split record_count into n_chunks sizes differing by at most one.

    The remainder goes to the leading chunks rather than piling onto the last,
    so no chunk is ever empty and none is disproportionately large.
    """
    base, remainder = divmod(record_count, n_chunks)
    return [base + (1 if index < remainder else 0) for index in range(n_chunks)]


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)

    input_vcf = sys.argv[1]
    requested_chunks = int(sys.argv[2])
    out_prefix = sys.argv[3]

    if requested_chunks < 1:
        sys.exit(f"n_chunks must be >= 1, got {requested_chunks}")

    header_lines, record_count = scan_vcf(input_vcf)

    if record_count == 0:
        sys.exit(f"No variant records found in {input_vcf}")
    if not any(line.startswith("#CHROM") for line in header_lines):
        sys.exit(f"No #CHROM line found in header of {input_vcf}")

    # Capping at record_count keeps every emitted chunk non-empty, which
    # cellsnp-lite may otherwise reject as a -T file.
    n_chunks = min(requested_chunks, record_count)
    sizes = chunk_sizes(record_count, n_chunks)

    # Second pass: walk the VCF in order, closing each chunk once it has taken
    # its quota. Only one output handle is open at a time, since the chunks are
    # contiguous.
    chunk_index = 0
    records_in_chunk = 0
    chunk_paths = []

    def open_chunk(index):
        path = f"{out_prefix}{index:03d}.vcf"
        handle = open(path, "wt")
        handle.writelines(header_lines)
        chunk_paths.append(path)
        return handle

    handle_out = open_chunk(chunk_index)

    with open_maybe_gz(input_vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            if records_in_chunk == sizes[chunk_index]:
                handle_out.close()
                chunk_index += 1
                handle_out = open_chunk(chunk_index)
                records_in_chunk = 0

            handle_out.write(line)
            records_in_chunk += 1

    handle_out.close()

    print(
        f"Split {record_count} SNP records into {len(chunk_paths)} contiguous "
        f"chunks of {sizes[-1]}-{sizes[0]} records (requested "
        f"{requested_chunks})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
