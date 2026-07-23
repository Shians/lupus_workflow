#!/usr/bin/env python3
"""Merge per-chunk cellsnp-lite output directories into one cellSNP directory.

Sharding the target SNP VCF and running cellsnp-lite independently per chunk is
exact: cellsnp-lite's ``--minMAF`` / ``--minCOUNT`` filters are per-SNP and use
only that SNP's own counts, and every chunk uses the *same* barcode file, so the
cell columns are identical and identically ordered across chunks. Merging is
therefore a concatenation along the SNP (row) axis.

The invariant that makes this correct: matrix row *i* corresponds to record *i*
of ``cellSNP.base.vcf`` and record *i* of ``cellSNP.cells.vcf``. vireo's
``read_cellSNP`` pairs matrix rows and VCF records by position only (it does not
sort), so all files must be concatenated in the *same* chunk order. Chunks are
ordered by directory basename, which is zero-padded with the chunk index, so
lexical order equals chunk order regardless of the order they arrive in.

Output filenames match what vireo's ``read_cellSNP`` expects, and the
compression of each output matches the corresponding per-chunk file (whatever
``cellsnp-lite --gzip`` produced), so the merged directory is byte-format
identical to a single monolithic cellsnp-lite run.

Usage:
    merge_cellsnp.py <out_dir> <chunk_dir> [<chunk_dir> ...]
"""

import gzip
import os
import shutil
import sys

# cellsnp-lite may or may not append ``.gz`` depending on --gzip and version;
# candidates are listed most-compressed-first and the present one is detected.
BASE_VCF = ["cellSNP.base.vcf.gz", "cellSNP.base.vcf"]
CELLS_VCF = ["cellSNP.cells.vcf.gz", "cellSNP.cells.vcf"]
MTX_FILES = {
    "AD": ["cellSNP.tag.AD.mtx.gz", "cellSNP.tag.AD.mtx"],
    "DP": ["cellSNP.tag.DP.mtx.gz", "cellSNP.tag.DP.mtx"],
    "OTH": ["cellSNP.tag.OTH.mtx.gz", "cellSNP.tag.OTH.mtx"],
}
SAMPLES_TSV = "cellSNP.samples.tsv"


def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def resolve(chunk_dir, candidates):
    """Return the path of the first existing candidate filename in chunk_dir."""
    for name in candidates:
        path = os.path.join(chunk_dir, name)
        if os.path.exists(path):
            return path
    raise FileNotFoundError(
        f"None of {candidates} found in chunk directory {chunk_dir}"
    )


def basename_of(chunk_dir, candidates):
    """Return the resolved output basename (preserving detected compression)."""
    return os.path.basename(resolve(chunk_dir, candidates))


def merge_samples(chunk_dirs, out_dir):
    """Copy samples.tsv from chunk 0; assert byte-identical across all chunks."""
    reference = os.path.join(chunk_dirs[0], SAMPLES_TSV)
    with open(reference, "rb") as handle:
        reference_bytes = handle.read()

    for chunk_dir in chunk_dirs[1:]:
        path = os.path.join(chunk_dir, SAMPLES_TSV)
        with open(path, "rb") as handle:
            if handle.read() != reference_bytes:
                sys.exit(
                    f"{SAMPLES_TSV} differs between {chunk_dirs[0]} and "
                    f"{chunk_dir}; chunks must share an identical barcode set"
                )

    shutil.copyfile(reference, os.path.join(out_dir, SAMPLES_TSV))
    # Number of barcodes = number of matrix columns; used to cross-check dims.
    return sum(1 for line in reference_bytes.splitlines() if line.strip())


def read_mtx_dims(path):
    """Return (header_lines, (nrows, ncols, nnz)) without loading the body."""
    header_lines = []
    with open_maybe_gz(path) as handle:
        for line in handle:
            if line.startswith("%"):
                header_lines.append(line)
                continue
            if not line.strip():
                continue
            parts = line.split()
            dims = (int(parts[0]), int(parts[1]), int(parts[2]))
            return header_lines, dims
    raise ValueError(f"No MatrixMarket dimension line found in {path}")


def merge_mtx(chunk_dirs, candidates, out_dir):
    """Concatenate one .mtx tag along the SNP (row) axis, re-offsetting rows."""
    out_name = basename_of(chunk_dirs[0], candidates)
    out_path = os.path.join(out_dir, out_name)
    open_out = gzip.open if out_name.endswith(".gz") else open

    # Pass 1: gather per-chunk dimensions; validate constant column count.
    header_lines = None
    ncols = None
    total_rows = 0
    total_nnz = 0
    per_chunk = []  # (resolved_path, nrows) in chunk order
    for chunk_dir in chunk_dirs:
        path = resolve(chunk_dir, candidates)
        chunk_header, (nrows, chunk_ncols, nnz) = read_mtx_dims(path)
        if header_lines is None:
            header_lines = chunk_header
        if ncols is None:
            ncols = chunk_ncols
        elif chunk_ncols != ncols:
            sys.exit(
                f"{out_name}: column count mismatch ({chunk_ncols} in "
                f"{chunk_dir} vs {ncols} in {chunk_dirs[0]}); barcode sets differ"
            )
        per_chunk.append((path, nrows))
        total_rows += nrows
        total_nnz += nnz

    # Pass 2: stream entries, adding the cumulative row offset for each chunk.
    with open_out(out_path, "wt") as out:
        out.writelines(header_lines)
        out.write(f"{total_rows} {ncols} {total_nnz}\n")

        row_offset = 0
        for path, nrows in per_chunk:
            seen_dims = False
            with open_maybe_gz(path) as handle:
                for line in handle:
                    if line.startswith("%"):
                        continue
                    if not line.strip():
                        continue
                    if not seen_dims:
                        seen_dims = True  # skip this chunk's dimension line
                        continue
                    parts = line.split()
                    row = int(parts[0]) + row_offset
                    col = parts[1]
                    value = parts[2]
                    out.write(f"{row} {col} {value}\n")
            row_offset += nrows

    return out_name, total_rows, ncols, total_nnz


def merge_vcf(chunk_dirs, candidates, out_dir):
    """Keep header from chunk 0; concatenate record lines in chunk order."""
    out_name = basename_of(chunk_dirs[0], candidates)
    out_path = os.path.join(out_dir, out_name)
    open_out = gzip.open if out_name.endswith(".gz") else open

    record_count = 0
    with open_out(out_path, "wt") as out:
        # Header (## meta + #CHROM) from chunk 0 only.
        with open_maybe_gz(resolve(chunk_dirs[0], candidates)) as handle:
            for line in handle:
                if line.startswith("#"):
                    out.write(line)

        # Record (non-#) lines from every chunk, in chunk order.
        for chunk_dir in chunk_dirs:
            with open_maybe_gz(resolve(chunk_dir, candidates)) as handle:
                for line in handle:
                    if not line.startswith("#"):
                        out.write(line)
                        record_count += 1

    return out_name, record_count


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)

    out_dir = sys.argv[1]
    chunk_dirs = sys.argv[2:]

    # Fix a single chunk order (by basename) and use it identically everywhere.
    chunk_dirs = sorted(chunk_dirs, key=os.path.basename)
    os.makedirs(out_dir, exist_ok=True)

    n_cells = merge_samples(chunk_dirs, out_dir)

    base_name, base_records = merge_vcf(chunk_dirs, BASE_VCF, out_dir)
    cells_name, cells_records = merge_vcf(chunk_dirs, CELLS_VCF, out_dir)
    if base_records != cells_records:
        sys.exit(
            f"Record count mismatch: {base_name} has {base_records} but "
            f"{cells_name} has {cells_records}"
        )

    for tag, candidates in MTX_FILES.items():
        mtx_name, mtx_rows, mtx_cols, mtx_nnz = merge_mtx(
            chunk_dirs, candidates, out_dir
        )
        if mtx_rows != base_records:
            sys.exit(
                f"Row/record mismatch: {mtx_name} has {mtx_rows} rows but "
                f"{base_name} has {base_records} records"
            )
        if mtx_cols != n_cells:
            sys.exit(
                f"Column/barcode mismatch: {mtx_name} has {mtx_cols} columns "
                f"but {SAMPLES_TSV} lists {n_cells} barcodes"
            )
        print(
            f"  {mtx_name}: {mtx_rows} SNPs x {mtx_cols} cells, {mtx_nnz} nnz",
            file=sys.stderr,
        )

    print(
        f"Merged {len(chunk_dirs)} chunks -> {out_dir}: "
        f"{base_records} SNPs x {n_cells} cells",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
