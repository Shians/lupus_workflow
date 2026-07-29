#!/usr/bin/env python3
"""Combine the per-sample depth tables from computeTargetDepth into one
kept-site list plus a per-site cost weight for splitCellSNPTargets.

Each input is a mosdepth ``regions.bed.gz`` over the candidate sites, gzipped:

    chrom  start  end  depth

``depth`` is CIGAR-aware coverage -- reads actually contributing a base at the
site, which is what can become an allele count. See the commentary above
computeTargetDepth for why that is the right measurement here and why
``--fast-mode`` is not.

Two aggregations, deliberately different:

* **keep** if ANY sample clears the floor. A site one sample covers is worth
  genotyping even if the others do not, because every sample is genotyped
  against one shared target list.
* **cost** is the MAX across samples, not the mean. Every sample pays its own
  pileup at every kept site, so the shard containing a site must be sized for
  the worst sample, not the average one.

Cost reuses the same depth rather than measuring pileup span separately. It is
a proxy: htslib's queue holds intron-spanning reads that contribute no base, so
depth understates the true queue. Among sites that survive the floor the two
track each other closely, and the weight only balances shards -- it changes no
result. See split_cellsnp_targets.py.

All inputs derive from the same candidate_sites.bed and so share a row order;
that is verified rather than assumed, and the outputs preserve it. Preserving
it is what lets splitCellSNPTargets stream the cost file alongside the filtered
VCF without holding either in memory.

``min_sites`` is a tripwire, not a quality bar. Demultiplexing degrades
gradually as sites are lost, so a run that keeps far too few does not fail --
it produces a Vireo assignment that looks ordinary and is quietly unreliable.
The floor turns that into an error. A healthy genome-wide run should clear it
by orders of magnitude; tripping it means something upstream is wrong (shallow
libraries, a reference or contig mismatch, a target list that does not match
the data) rather than that the floor needs lowering.

Usage:
    aggregate_target_depth.py <min_depth> <min_sites> <kept.bed> <cost.tsv> \
        <depth.bed.gz>...
"""

import gzip
import sys


def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse(line, path, line_number):
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 4:
        sys.exit(
            f"{path} line {line_number}: expected 4 tab-separated fields "
            f"(chrom, start, end, depth), got {len(fields)}"
        )
    chrom, start, end, depth = fields
    try:
        # mosdepth reports a regional mean, so this is a float even for the
        # 1 bp regions used here.
        return chrom, start, end, float(depth)
    except ValueError:
        sys.exit(f"{path} line {line_number}: depth column is not numeric")


def main():
    if len(sys.argv) < 6:
        sys.exit(__doc__)

    min_depth = float(sys.argv[1])
    min_sites = int(sys.argv[2])
    kept_path = sys.argv[3]
    cost_path = sys.argv[4]
    depth_paths = sys.argv[5:]

    handles = [open_maybe_gz(path) for path in depth_paths]
    kept_count = 0
    site_count = 0

    with open(kept_path, "wt") as kept, open(cost_path, "wt") as cost:
        for line_number, rows in enumerate(zip(*handles), start=1):
            site_count += 1
            chrom = start = end = None
            keep = False
            max_depth = 0.0

            for path, line in zip(depth_paths, rows):
                row_chrom, row_start, row_end, depth = parse(
                    line, path, line_number
                )
                if chrom is None:
                    chrom, start, end = row_chrom, row_start, row_end
                elif (row_chrom, row_start) != (chrom, start):
                    sys.exit(
                        f"line {line_number}: {depth_paths[0]} has "
                        f"{chrom}:{start} but {path} has {row_chrom}:{row_start}; "
                        f"per-sample depth tables must share a row order"
                    )
                if depth >= min_depth:
                    keep = True
                max_depth = max(max_depth, depth)

            if keep:
                kept_count += 1
                kept.write(f"{chrom}\t{start}\t{end}\n")
                # Floor of 1 so that a kept site never has zero weight: cost
                # balancing still has to spread the sites themselves, not only
                # the depth they carry.
                cost.write(f"{chrom}\t{end}\t{max(1, round(max_depth))}\n")

        # zip() stops at the shortest input, so a truncated file would silently
        # shorten the output rather than fail. Catch that explicitly.
        for path, handle in zip(depth_paths, handles):
            if handle.readline():
                sys.exit(
                    f"{path} has more rows than the shortest input, which "
                    f"stopped at line {site_count}; per-sample depth tables "
                    f"must cover the identical site list"
                )

    for handle in handles:
        handle.close()

    if kept_count < min_sites:
        sys.exit(
            f"Only {kept_count} of {site_count} target sites reached depth "
            f"{min_depth} in any sample, below the {min_sites} required to "
            f"genotype reliably. Demultiplexing on this few sites would not "
            f"fail outright, it would just be untrustworthy, so this stops "
            f"here. Check library depth, and that the target VCF's contigs "
            f"and coordinates match the aligned reference, before lowering "
            f"min_target_depth or min_target_sites."
        )

    print(
        f"Kept {kept_count} of {site_count} target sites at depth "
        f">= {min_depth} in at least one sample",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
