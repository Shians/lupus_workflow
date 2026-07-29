#!/usr/bin/env python3
"""Split a cellSNP target SNP VCF into N contiguous, cost-balanced chunks.

Each chunk carries the full VCF header (all ``##`` meta lines plus the
``#CHROM`` line) so that cellsnp-lite accepts it as a ``-T`` target file.

Two properties matter, and the obvious split satisfies only one of them:

**Data locality.** Chunks are *contiguous* runs of the coordinate-sorted target
list, so each shard reads one genomic interval of the BAM. The previous
round-robin split gave every shard sites on every chromosome, which meant every
shard traversed the whole BAM -- 64 shards over a 276 GB BAM is ~17.7 TB of
reads to genotype what one pass covers in 276 GB. Contiguity removes that
amplification entirely.

**Balance.** Contiguous slicing with equal record counts is badly unbalanced,
which is why round-robin was chosen originally: cost per site spans orders of
magnitude, and equal counts hand one shard a cluster of pathologically deep
sites. So chunks are balanced on summed *cost* instead, where cost is the
per-site read depth measured by computeTargetDepth and aggregated across
samples by aggregate_target_depth.py. That recovers the balance round-robin was
providing without giving up locality.

Depth is a proxy for what a site actually costs, which is the size of htslib's
pileup queue there -- and the queue also holds intron-spanning reads that
contribute no base, so depth understates it. Measuring the queue directly means
a second, span-based mosdepth pass for a number that balances shards and
changes no result. Not worth it: the two diverge sharply only under introns,
and those sites no longer survive the depth floor. A mis-weighted site just
makes its shard run long, and cellsnp_max_pileup bounds the tail regardless.

Note the consequence for threading: a contiguous chunk usually spans a single
chromosome, and cellsnp-lite parallelises by chromosome, so a shard cannot use
many threads. That is intended -- shards are meant to be narrow and numerous,
and one pileup queue per shard rather than one per thread is a large part of
the memory saving. Size the chunk tier's cpus accordingly.

The cost file is optional. When absent or empty every record weighs 1, which
degrades to a contiguous equal-count split -- still local, just less balanced.

Usage:
    split_cellsnp_targets.py <input.vcf[.gz]> <n_chunks> <out_prefix> [cost.tsv]

Emits files ``<out_prefix>000.vcf``, ``<out_prefix>001.vcf``, ... one per
non-empty chunk. Output is plain (uncompressed) VCF, which cellsnp-lite reads
directly with no index required.
"""

import gzip
import os
import sys


def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_costs(cost_path):
    """Yield ``(chrom, pos, cost)`` from the aggregate_target_depth.py table."""
    with open_maybe_gz(cost_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3:
                sys.exit(
                    f"{cost_path} line {line_number}: expected 3 tab-separated "
                    f"fields (chrom, pos, cost), got {len(fields)}"
                )
            try:
                yield fields[0], fields[1], float(fields[2])
            except ValueError:
                sys.exit(f"{cost_path} line {line_number}: cost is not numeric")


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


def main():
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)

    input_vcf = sys.argv[1]
    requested_chunks = int(sys.argv[2])
    out_prefix = sys.argv[3]
    cost_path = sys.argv[4] if len(sys.argv) == 5 else None

    if requested_chunks < 1:
        sys.exit(f"n_chunks must be >= 1, got {requested_chunks}")

    if cost_path is not None and (
        not os.path.exists(cost_path) or os.path.getsize(cost_path) == 0
    ):
        cost_path = None

    header_lines, record_count = scan_vcf(input_vcf)

    if record_count == 0:
        sys.exit(f"No variant records found in {input_vcf}")
    if not any(line.startswith("#CHROM") for line in header_lines):
        sys.exit(f"No #CHROM line found in header of {input_vcf}")

    n_chunks = min(requested_chunks, record_count)

    # First pass over the costs, to know the total before deciding where the
    # boundaries fall. Streamed twice rather than held in memory: the target
    # list runs to tens of millions of sites.
    if cost_path is None:
        total_cost = float(record_count)
    else:
        total_cost = 0.0
        cost_records = 0
        for _chrom, _pos, cost in read_costs(cost_path):
            total_cost += cost
            cost_records += 1
        if cost_records != record_count:
            sys.exit(
                f"{cost_path} has {cost_records} rows but {input_vcf} has "
                f"{record_count} records; the cost table must describe exactly "
                f"the sites being split"
            )

    # Second pass: walk the VCF and the costs in lockstep, closing the current
    # chunk once it has taken its share of the cost still unallocated.
    #
    # The target is recomputed from what remains after every cut rather than
    # fixed at total/n_chunks up front, because a single site can cost more than
    # a whole chunk's share -- the deepest here run to millions while the median
    # is 13. Against fixed thresholds one such site overshoots several
    # boundaries at once, and since a cut advances the chunk index only one
    # step, the overshoot drains off as a run of one-site chunks. Rebalancing
    # absorbs the overshoot into the remaining budget instead, so an expensive
    # site simply takes a chunk to itself and the rest stay even.
    costs = read_costs(cost_path) if cost_path is not None else None
    chunk_index = 0
    chunk_cost = 0.0
    records_in_chunk = 0
    records_written = 0
    remaining_cost = total_cost
    remaining_chunks = n_chunks
    target = remaining_cost / remaining_chunks
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

            if costs is None:
                cost = 1.0
            else:
                chrom, pos, cost = next(costs)
                vcf_chrom, vcf_pos = line.split("\t", 2)[:2]
                if (chrom, pos) != (vcf_chrom, vcf_pos):
                    sys.exit(
                        f"cost table is at {chrom}:{pos} but {input_vcf} is at "
                        f"{vcf_chrom}:{vcf_pos}; the two must be in the same "
                        f"coordinate order over the same sites"
                    )

            # Decide before writing, so a boundary starts the new chunk with
            # this record rather than stranding it at the end of the old one.
            records_left = record_count - records_written
            if remaining_chunks > 1 and records_in_chunk > 0:
                # Cut when this chunk has taken its share, or when the records
                # left are only just enough to give every remaining chunk one.
                # The second clause is what guarantees no chunk comes out
                # empty, which cellsnp-lite may reject as a -T file.
                if chunk_cost >= target or records_left <= remaining_chunks:
                    handle_out.close()
                    remaining_cost -= chunk_cost
                    remaining_chunks -= 1
                    target = remaining_cost / remaining_chunks
                    chunk_index += 1
                    handle_out = open_chunk(chunk_index)
                    chunk_cost = 0.0
                    records_in_chunk = 0

            handle_out.write(line)
            chunk_cost += cost
            records_in_chunk += 1
            records_written += 1

    handle_out.close()

    if costs is not None and next(costs, None) is not None:
        sys.exit(f"{cost_path} has more rows than {input_vcf} has records")

    print(
        f"Split {record_count} SNP records into {len(chunk_paths)} contiguous "
        f"chunks (requested {requested_chunks}), balanced by "
        f"{'depth cost' if cost_path else 'record count'}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
