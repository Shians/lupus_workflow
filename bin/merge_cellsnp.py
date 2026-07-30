#!/usr/bin/env python3
"""Merge per-chunk cellsnp-lite output directories into one cellSNP directory.

Sharding the target SNP VCF and running cellsnp-lite independently per chunk is
exact: cellsnp-lite's ``--minMAF`` / ``--minCOUNT`` filters are per-SNP and use
only that SNP's own counts, and every chunk uses the *same* barcode file, so the
cell columns are identical and identically ordered across chunks. Merging
therefore only has to get the SNP (row) axis right; the cell columns need no
attention at all.

``split_cellsnp_targets.py`` emits *contiguous* ranges of the coordinate-sorted
target VCF: chunk 000 holds the first N records, chunk 001 the next N, and so
on. cellsnp-lite emits sites in target-file order. So the merged order is plain
concatenation of the chunks in index order, and a record's merged row is its
local row plus a per-chunk offset -- an integer add rather than a permutation.
No heap, no coordinate comparison, no row map.

That makes this script *depend* on the split being contiguous, which the k-way
merge it replaces did not. Reverting the splitter to round-robin would silently
produce a mis-ordered directory, so the coupling is checked rather than trusted:
``cellSNP.base.vcf`` is verified non-decreasing within each chunk, each chunk's
last coordinate is verified to precede the next chunk's first, and every ``.mtx``
is verified row-major. A tripped check fails loudly.

The invariant that makes the result usable: matrix row *i* must correspond to
record *i* of ``cellSNP.base.vcf``. vireo's ``read_cellSNP`` pairs matrix rows
with VCF records by position only (it does not sort), so one row numbering is
applied to the base VCF and to all three matrices alike -- and to
``cellSNP.cells.vcf`` too, on the rare runs that have one.

``cellSNP.cells.vcf`` is *optional* here and normally absent. It exists only
under ``cellsnp-lite --genotype``, which the pipeline no longer passes because
no consumer reads the file: vireo's ``read_cellSNP`` and
``snplet::import_cellsnp`` both open only ``base.vcf``, the matrices and
``samples.tsv``. When it is present it is copied as raw blocks rather than line
by line -- its records carry one field per cell, so a line-oriented pass would
build ~200 kB Python strings per record for no purpose. Row ordering is settled
by ``cellSNP.base.vcf``, which is small and parsed in full.

**On compression.** In the default configuration this is a minor cost: only
``cellSNP.base.vcf.gz`` is gzipped, since cellsnp-lite writes the ``.mtx`` files
as plain text, so the matrix entry loop dominates instead. The numbers below are
kept because they are why ``--genotype`` was dropped, and because they apply
again to anyone who restores it.

With ``--genotype``, ``cellSNP.cells.vcf.gz`` inflates 167x at 41,259 cells
(33 MB per shard on disk, 5.5 GB of text), so merging it moved ~420 GB through
Python's ``gzip`` twice. That cost is invisible in a Nextflow trace:
``rchar``/``wchar`` count compressed bytes at the syscall boundary, so such a
file contributes its *compressed* size to the counters and its *uncompressed*
size to the runtime. Measured on real data (844,477 SNPs x 41,259 cells),
decompression ran at 460 MB/s against 96 MB/s for the stdlib's default level 9
-- ~64 min of a 67 min task, on a file nothing opens.

Hence two choices that survive regardless. Level 6 rather than 9 (182 MB/s vs
96 MB/s on that data, for a few percent of size), and ``--threads`` handed to an
external compressor when one is on PATH. ``bgzip`` is preferred over ``pigz``:
cellsnp-lite emits BGZF, and recompressing with plain gzip silently costs the
block structure that makes the file indexable. All three paths produce valid
gzip, so correctness does not depend on which is found.

Output filenames match what vireo's ``read_cellSNP`` expects, and the
compression of each output matches the corresponding per-chunk file (whatever
``cellsnp-lite --gzip`` produced -- in practice gzipped VCFs and plain ``.mtx``).
"""

import argparse
import gzip
import os
import shutil
import subprocess
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

COPY_BLOCK = 8 << 20


def open_read(path):
    """Open a possibly-gzipped file for binary reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def compressor_command(level, threads):
    """Return the argv of an external stdin->stdout gzip compressor, or None.

    bgzip first, and not only for its threading: cellsnp-lite writes BGZF, so
    recompressing with anything else silently downgrades the merged file to
    plain gzip and loses the block structure that makes it indexable. pigz is
    the threaded fallback, and the stdlib is the correct-but-serial floor.
    """
    bgzip = shutil.which("bgzip")
    if bgzip is not None:
        return [bgzip, "-l", str(level), "-@", str(threads), "-c"]

    pigz = shutil.which("pigz")
    if pigz is not None:
        return [pigz, f"-{level}", "-p", str(threads)]

    return None


class OutputStream:
    """Binary writer that gzips through an external compressor when available."""

    def __init__(self, path, level, threads):
        self.path = path
        self.process = None
        self.raw_file = None
        self.command = None

        if not path.endswith(".gz"):
            self.handle = open(path, "wb")
            return

        command = compressor_command(level, threads)
        if command is None:
            self.handle = gzip.open(path, "wb", compresslevel=level)
            return

        self.command = command
        self.raw_file = open(path, "wb")
        self.process = subprocess.Popen(
            command, stdin=subprocess.PIPE, stdout=self.raw_file
        )
        self.handle = self.process.stdin

    def __enter__(self):
        return self.handle

    def __exit__(self, *exc_info):
        self.handle.close()
        if self.process is not None:
            returncode = self.process.wait()
            self.raw_file.close()
            if returncode != 0 and exc_info[0] is None:
                sys.exit(
                    f"{os.path.basename(self.command[0])} failed with status "
                    f"{returncode} writing {self.path}"
                )
        return False


def find(chunk_dir, candidates):
    """Return the path of the first existing candidate filename, or None."""
    for name in candidates:
        path = os.path.join(chunk_dir, name)
        if os.path.exists(path):
            return path
    return None


def resolve(chunk_dir, candidates):
    """Return the path of the first existing candidate filename in chunk_dir."""
    path = find(chunk_dir, candidates)
    if path is None:
        raise FileNotFoundError(
            f"None of {candidates} found in chunk directory {chunk_dir}"
        )
    return path


def basename_of(chunk_dir, candidates):
    """Return the resolved output basename (preserving detected compression)."""
    return os.path.basename(resolve(chunk_dir, candidates))


def read_contig_ranks(path):
    """Map contig name -> rank from the sorted target VCF's contig order.

    Taken from the order contigs actually appear in the sorted VCF rather than
    from its ``##contig`` header lines, which need not be present, complete, or
    in agreement with the records. Keyed by bytes, since the VCFs are read in
    binary.
    """
    ranks = {}
    with open(path, "rb") as handle:
        for line in handle:
            name = line.strip()
            if name and name not in ranks:
                ranks[name] = len(ranks)
    if not ranks:
        sys.exit(f"No contig names found in {path}")
    return ranks


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


def split_header(handle):
    """Consume the leading ``#`` lines, returning them and the first record."""
    header = []
    while True:
        line = handle.readline()
        if not line:
            return header, None
        if line.startswith(b"#"):
            header.append(line)
        else:
            return header, line


def record_key(line, contig_ranks, chunk_dir):
    """Extract ``(contig_rank, pos)`` from a VCF record without splitting it."""
    first_tab = line.find(b"\t")
    second_tab = line.find(b"\t", first_tab + 1)
    if first_tab < 0 or second_tab < 0:
        sys.exit(f"{chunk_dir}: malformed VCF record: {line[:80]!r}")

    contig = line[:first_tab]
    rank = contig_ranks.get(contig)
    if rank is None:
        sys.exit(
            f"{chunk_dir}: contig {contig.decode(errors='replace')!r} is absent "
            f"from the target VCF contig order; this chunk was not produced "
            f"from that VCF"
        )
    return rank, int(line[first_tab + 1:second_tab])


def concat_base_vcf(chunk_dirs, contig_ranks, out_dir, level, threads):
    """Concatenate ``cellSNP.base.vcf``, verifying the chunks are in order.

    Returns ``(record_counts, out_name)`` where ``record_counts[i]`` is chunk
    *i*'s record count. Their running sum gives the offset each chunk's matrix
    rows are shifted by.

    This is the only file whose records are parsed. It carries no per-cell
    fields, so the full coordinate check is cheap here, and every other file in
    the directory is row-parallel to it and inherits the result.
    """
    out_name = basename_of(chunk_dirs[0], BASE_VCF)
    record_counts = []
    previous_chunk_last = None
    previous_chunk_dir = None

    with OutputStream(os.path.join(out_dir, out_name), level, threads) as out:
        for chunk_dir in chunk_dirs:
            with open_read(resolve(chunk_dir, BASE_VCF)) as handle:
                header, line = split_header(handle)
                # Headers come from chunk 0; every chunk saw the same BAM and
                # barcodes, so the others' headers are redundant.
                if not record_counts:
                    out.writelines(header)

                count = 0
                previous_key = None
                while line:
                    key = record_key(line, contig_ranks, chunk_dir)

                    if previous_key is not None and key < previous_key:
                        sys.exit(
                            f"{chunk_dir}: records are not in coordinate order "
                            f"({previous_key} then {key}). cellsnp-lite did not "
                            f"emit sites in target-file order, so concatenating "
                            f"the chunks cannot reproduce a sorted directory"
                        )
                    if count == 0 and previous_chunk_last is not None \
                            and key < previous_chunk_last:
                        sys.exit(
                            f"{chunk_dir} starts at {key} but {previous_chunk_dir} "
                            f"ended at {previous_chunk_last}; the chunks are not "
                            f"contiguous ranges of the sorted target VCF, so "
                            f"concatenating them would mis-order the output. "
                            f"Did split_cellsnp_targets.py stop splitting into "
                            f"contiguous ranges?"
                        )

                    out.write(line)
                    if not line.endswith(b"\n"):
                        out.write(b"\n")
                    previous_key = key
                    count += 1
                    line = handle.readline()

            if count == 0:
                sys.exit(f"{chunk_dir}: {out_name} holds no records")

            record_counts.append(count)
            previous_chunk_last = previous_key
            previous_chunk_dir = chunk_dir

    return record_counts, out_name


def concat_cells_vcf(chunk_dirs, out_dir, record_counts, level, threads):
    """Concatenate ``cellSNP.cells.vcf`` as raw blocks, counting records.

    Returns ``None`` when the chunks carry no cells VCF, which is the normal
    case: the file exists only under ``cellsnp-lite --genotype`` and nothing
    downstream reads it. Absence must be unanimous -- a mix means the chunks
    came from different cellsnp-lite invocations, which puts every other
    cross-chunk assumption here in doubt.

    Copied block-wise rather than line-wise: each record carries one field per
    cell, so materialising it as a Python string costs far more than the newline
    count this needs. Ordering is not re-derived here -- it is settled by
    ``cellSNP.base.vcf``, and the per-chunk record counts are checked against it.
    """
    present = [find(d, CELLS_VCF) is not None for d in chunk_dirs]
    if not any(present):
        return None
    if not all(present):
        missing = [d for d, ok in zip(chunk_dirs, present) if not ok]
        sys.exit(
            f"{CELLS_VCF[0]} is present in some chunks but absent from "
            f"{len(missing)} of {len(chunk_dirs)} (e.g. {missing[0]}); the "
            f"chunks were not produced by the same cellsnp-lite invocation"
        )

    out_name = basename_of(chunk_dirs[0], CELLS_VCF)

    with OutputStream(os.path.join(out_dir, out_name), level, threads) as out:
        for chunk_index, chunk_dir in enumerate(chunk_dirs):
            with open_read(resolve(chunk_dir, CELLS_VCF)) as handle:
                header, line = split_header(handle)
                if chunk_index == 0:
                    out.writelines(header)

                count = 0
                last_byte = b""
                if line:
                    out.write(line)
                    count = 1
                    last_byte = line[-1:]

                while True:
                    block = handle.read(COPY_BLOCK)
                    if not block:
                        break
                    out.write(block)
                    count += block.count(b"\n")
                    last_byte = block[-1:]

                # A chunk whose final record lacks its newline would otherwise
                # be spliced onto the next chunk's first record.
                if last_byte and last_byte != b"\n":
                    out.write(b"\n")
                    count += 1

            if count != record_counts[chunk_index]:
                sys.exit(
                    f"{chunk_dir}: {out_name} holds {count} records but "
                    f"{basename_of(chunk_dir, BASE_VCF)} holds "
                    f"{record_counts[chunk_index]}; the two must be row-for-row "
                    f"parallel"
                )

    return out_name


def read_mtx_dims(path):
    """Return ``(header_lines, (rows, cols, nnz))`` without reading the entries."""
    header = []
    with open_read(path) as handle:
        while True:
            line = handle.readline()
            if not line:
                sys.exit(f"No MatrixMarket dimension line found in {path}")
            if line.startswith(b"%"):
                header.append(line)
                continue
            if not line.strip():
                continue
            parts = line.split()
            return header, (int(parts[0]), int(parts[1]), int(parts[2]))


def concat_mtx(chunk_dirs, candidates, out_dir, record_counts, level, threads):
    """Concatenate one .mtx tag, shifting each chunk's row indices by its offset.

    Row-major order is preserved for free: chunk *i*'s rows all sort below chunk
    *i+1*'s once the offsets are applied, and each chunk is already row-major.
    """
    out_name = basename_of(chunk_dirs[0], candidates)
    total_rows = sum(record_counts)

    # The output's dimension line needs the totals before any entry can be
    # written. Only each chunk's header is read here, not its entries.
    header_lines = None
    ncols = None
    total_nnz = 0
    for chunk_index, chunk_dir in enumerate(chunk_dirs):
        header, (nrows, chunk_ncols, nnz) = read_mtx_dims(
            resolve(chunk_dir, candidates)
        )
        if header_lines is None:
            header_lines = header
        if ncols is None:
            ncols = chunk_ncols
        elif chunk_ncols != ncols:
            sys.exit(
                f"{out_name}: column count mismatch ({chunk_ncols} in "
                f"{chunk_dir} vs {ncols} in {chunk_dirs[0]}); barcode sets differ"
            )
        if nrows != record_counts[chunk_index]:
            sys.exit(
                f"{out_name}: {chunk_dir} declares {nrows} rows but its VCFs "
                f"hold {record_counts[chunk_index]} records"
            )
        total_nnz += nnz

    written_nnz = 0
    offset = 0
    with OutputStream(os.path.join(out_dir, out_name), level, threads) as out:
        out.writelines(header_lines)
        out.write(f"{total_rows} {ncols} {total_nnz}\n".encode())

        for chunk_index, chunk_dir in enumerate(chunk_dirs):
            path = resolve(chunk_dir, candidates)
            nrows = record_counts[chunk_index]
            previous_row = 0

            with open_read(path) as handle:
                # Skip the header and dimension line already accounted for.
                while True:
                    line = handle.readline()
                    if not line:
                        break
                    if line.startswith(b"%") or not line.strip():
                        continue
                    break

                while True:
                    line = handle.readline()
                    if not line:
                        break
                    if not line.strip():
                        continue

                    # Only the row index is rewritten; the column and value are
                    # passed through as the bytes they arrived as.
                    field_end = line.find(b" ")
                    if field_end < 0:
                        sys.exit(f"{path}: malformed .mtx entry: {line[:80]!r}")
                    row = int(line[:field_end])

                    if row < previous_row:
                        sys.exit(
                            f"{path}: entries are not row-major (row {row} "
                            f"appears after row {previous_row}); rows must "
                            f"ascend for the offset shift to preserve ordering"
                        )
                    if row < 1 or row > nrows:
                        sys.exit(
                            f"{path}: row index {row} is outside the declared "
                            f"1..{nrows} range"
                        )

                    out.write(str(row + offset).encode())
                    out.write(line[field_end:])
                    if not line.endswith(b"\n"):
                        out.write(b"\n")
                    previous_row = row
                    written_nnz += 1

            offset += nrows

    if written_nnz != total_nnz:
        sys.exit(
            f"{out_name}: wrote {written_nnz} entries but the chunk headers "
            f"declared {total_nnz}"
        )

    return out_name, total_rows, ncols, total_nnz


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Concatenate contiguous per-chunk cellsnp-lite output directories "
            "into one cellSNP directory."
        )
    )
    parser.add_argument("out_dir", help="directory to write the merged output to")
    parser.add_argument(
        "contig_order",
        help="contig names in the order the sorted target VCF used, one per line",
    )
    parser.add_argument(
        "chunk_dirs", nargs="+", help="per-chunk cellsnp-lite output directories"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help=(
            "threads for bgzip, or pigz if bgzip is absent; ignored when neither "
            "is on PATH and the stdlib does the work (default: %(default)s)"
        ),
    )
    parser.add_argument(
        "--compress-level",
        type=int,
        default=6,
        help=(
            "gzip level for the gzipped outputs. The stdlib default of 9 measured "
            "roughly half the throughput of 6 on real data, for a few percent of "
            "size (default: %(default)s)"
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Sorted because chunk index order *is* the coordinate order: chunk 000
    # holds the first records of the sorted target VCF. concat_base_vcf checks
    # that this holds rather than assuming it.
    chunk_dirs = sorted(args.chunk_dirs, key=os.path.basename)

    os.makedirs(args.out_dir, exist_ok=True)
    contig_ranks = read_contig_ranks(args.contig_order)

    n_cells = merge_samples(chunk_dirs, args.out_dir)

    record_counts, base_name = concat_base_vcf(
        chunk_dirs, contig_ranks, args.out_dir, args.compress_level, args.threads
    )
    total_records = sum(record_counts)
    if total_records == 0:
        sys.exit(f"No records found across {len(chunk_dirs)} chunk directories")

    cells_name = concat_cells_vcf(
        chunk_dirs, args.out_dir, record_counts, args.compress_level, args.threads
    )
    written = base_name if cells_name is None else f"{base_name} / {cells_name}"
    print(f"  {written}: {total_records} SNPs", file=sys.stderr)

    for candidates in MTX_FILES.values():
        mtx_name, mtx_rows, mtx_cols, mtx_nnz = concat_mtx(
            chunk_dirs,
            candidates,
            args.out_dir,
            record_counts,
            args.compress_level,
            args.threads,
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
        f"Merged {len(chunk_dirs)} chunks -> {args.out_dir}: "
        f"{total_records} SNPs x {n_cells} cells, coordinate-sorted",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
