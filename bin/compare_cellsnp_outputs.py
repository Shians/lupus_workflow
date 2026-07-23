#!/usr/bin/env python3
"""Validate that a sharded+merged cellSNP directory equals a single-run one.

Implements the Task 6 assertions: two cellsnp-lite output directories must
describe the same SNP set (keyed by CHROM,POS,REF,ALT, order may differ) and,
after aligning rows by that key, have identical AD/DP/OTH matrices and identical
samples.tsv. Exits non-zero (with a diagnostic) on the first mismatch.

Usage:
    compare_cellsnp_outputs.py <cellsnp_dir_a> <cellsnp_dir_b>

Run this on a single test chromosome (or small SNP subset) before trusting the
sharded path at full 100-way fan-out. After it passes, also run vireo on both
directories and confirm identical donor_ids.tsv assignments.
"""

import gzip
import os
import sys

BASE_VCF = ["cellSNP.base.vcf.gz", "cellSNP.base.vcf"]
SAMPLES_TSV = "cellSNP.samples.tsv"
MTX_FILES = {
    "AD": ["cellSNP.tag.AD.mtx.gz", "cellSNP.tag.AD.mtx"],
    "DP": ["cellSNP.tag.DP.mtx.gz", "cellSNP.tag.DP.mtx"],
    "OTH": ["cellSNP.tag.OTH.mtx.gz", "cellSNP.tag.OTH.mtx"],
}


def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def resolve(directory, candidates):
    for name in candidates:
        path = os.path.join(directory, name)
        if os.path.exists(path):
            return path
    sys.exit(f"None of {candidates} found in {directory}")


def read_snp_keys(directory):
    """Return list of (CHROM,POS,REF,ALT) keys in base.vcf record order."""
    keys = []
    with open_maybe_gz(resolve(directory, BASE_VCF)) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            keys.append((fields[0], fields[1], fields[3], fields[4]))
    return keys


def read_mtx(path):
    """Return {(row, col): value} (1-based indices) and (nrows, ncols)."""
    entries = {}
    dims = None
    with open_maybe_gz(path) as handle:
        for line in handle:
            if line.startswith("%") or not line.strip():
                continue
            parts = line.split()
            if dims is None:
                dims = (int(parts[0]), int(parts[1]))
                continue
            entries[(int(parts[0]), int(parts[1]))] = parts[2]
    return entries, dims


def compare_matrix(tag, candidates, dir_a, dir_b, keys_a, keys_b):
    entries_a, dims_a = read_mtx(resolve(dir_a, candidates))
    entries_b, dims_b = read_mtx(resolve(dir_b, candidates))

    if dims_a[1] != dims_b[1]:
        sys.exit(f"{tag}: column counts differ ({dims_a[1]} vs {dims_b[1]})")

    # Map each SNP key to its 1-based row index in each directory.
    row_of_a = {key: i + 1 for i, key in enumerate(keys_a)}
    row_of_b = {key: i + 1 for i, key in enumerate(keys_b)}

    # Re-key both matrices by (snp_key, col) so row-order differences wash out.
    def rekey(entries, keys):
        rekeyed = {}
        for (row, col), value in entries.items():
            rekeyed[(keys[row - 1], col)] = value
        return rekeyed

    rekeyed_a = rekey(entries_a, keys_a)
    rekeyed_b = rekey(entries_b, keys_b)

    if rekeyed_a.keys() != rekeyed_b.keys():
        only_a = len(set(rekeyed_a) - set(rekeyed_b))
        only_b = len(set(rekeyed_b) - set(rekeyed_a))
        sys.exit(
            f"{tag}: non-zero entry pattern differs "
            f"({only_a} only in A, {only_b} only in B)"
        )

    for cell_key in rekeyed_a:
        if rekeyed_a[cell_key] != rekeyed_b[cell_key]:
            sys.exit(
                f"{tag}: value differs at {cell_key}: "
                f"{rekeyed_a[cell_key]} vs {rekeyed_b[cell_key]}"
            )

    _ = (row_of_a, row_of_b)  # row maps kept for clarity/debugging
    print(f"  {tag}: {len(rekeyed_a)} non-zero entries match")


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)

    dir_a, dir_b = sys.argv[1], sys.argv[2]

    # samples.tsv identical
    with open(os.path.join(dir_a, SAMPLES_TSV), "rb") as handle:
        samples_a = handle.read()
    with open(os.path.join(dir_b, SAMPLES_TSV), "rb") as handle:
        samples_b = handle.read()
    if samples_a != samples_b:
        sys.exit(f"{SAMPLES_TSV} differs between {dir_a} and {dir_b}")
    print(f"{SAMPLES_TSV}: identical")

    keys_a = read_snp_keys(dir_a)
    keys_b = read_snp_keys(dir_b)
    if set(keys_a) != set(keys_b):
        only_a = len(set(keys_a) - set(keys_b))
        only_b = len(set(keys_b) - set(keys_a))
        sys.exit(f"SNP sets differ ({only_a} only in A, {only_b} only in B)")
    if len(keys_a) != len(set(keys_a)):
        sys.exit(f"{dir_a}: duplicate SNP keys in base.vcf")
    if len(keys_b) != len(set(keys_b)):
        sys.exit(f"{dir_b}: duplicate SNP keys in base.vcf")
    print(f"SNP set: {len(keys_a)} SNPs match by (CHROM,POS,REF,ALT)")

    for tag, candidates in MTX_FILES.items():
        compare_matrix(tag, candidates, dir_a, dir_b, keys_a, keys_b)

    print("PASS: sharded+merged output matches single-run output")


if __name__ == "__main__":
    main()
