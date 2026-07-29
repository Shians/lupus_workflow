#!/usr/bin/env python3
"""Reduce the per-sample depth tables from computeTargetDepth to one kept-site
list for filterTargetsByDepth.

Each input is a mosdepth ``regions.bed.gz`` over the candidate sites, gzipped:

    chrom  start  end  depth

``depth`` is CIGAR-aware coverage -- reads actually contributing a base at the
site, which is what can become an allele count. See the commentary above
computeTargetDepth for why that is the right measurement here and why
``--fast-mode`` is not.

A site is kept if ANY sample clears the floor. A site one sample covers is
worth genotyping even if the others do not, because every sample is genotyped
against one shared target list. Per-sample-aware filtering still happens
downstream, via the per-SNP minMAF/minCOUNT filters in cellsnpDonorFilters.

All inputs derive from the same candidate_sites.bed and so share a row order,
which lets this read them row-for-row in lockstep and hold nothing in memory.
That is verified rather than assumed.

``min_sites`` is a tripwire, not a quality bar. Demultiplexing degrades
gradually as sites are lost, so a run that keeps far too few does not fail --
it produces a Vireo assignment that looks ordinary and is quietly unreliable.
The floor turns that into an error. A healthy genome-wide run should clear it
by orders of magnitude; tripping it means something upstream is wrong (shallow
libraries, a reference or contig mismatch, a target list that does not match
the data) rather than that the floor needs lowering.

Usage:
    aggregate_target_depth.py <min_depth> <min_sites> <kept.bed> <depth.bed.gz>...
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
    if len(sys.argv) < 5:
        sys.exit(__doc__)

    min_depth = float(sys.argv[1])
    min_sites = int(sys.argv[2])
    kept_path = sys.argv[3]
    depth_paths = sys.argv[4:]

    handles = [open_maybe_gz(path) for path in depth_paths]
    kept_count = 0
    site_count = 0

    with open(kept_path, "wt") as kept:
        for line_number, rows in enumerate(zip(*handles), start=1):
            site_count += 1
            chrom = start = end = None
            keep = False

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

            if keep:
                kept_count += 1
                kept.write(f"{chrom}\t{start}\t{end}\n")

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
