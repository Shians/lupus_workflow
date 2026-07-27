#!/usr/bin/env python3
"""Merge per-chunk cellsnp-lite output directories into one cellSNP directory.

Sharding the target SNP VCF and running cellsnp-lite independently per chunk is
exact: cellsnp-lite's ``--minMAF`` / ``--minCOUNT`` filters are per-SNP and use
only that SNP's own counts, and every chunk uses the *same* barcode file, so the
cell columns are identical and identically ordered across chunks. Merging is
therefore a reordering along the SNP (row) axis.

The target VCF is coordinate-sorted before it is sharded (see
``sortSNPAnnotation``), and cellsnp-lite emits sites in target-file order, so
each chunk's output is itself coordinate-sorted: a round-robin subset of a
sorted list stays sorted. Merging is therefore a k-way merge over already
sorted streams, which reconstructs the order a single monolithic run would have
produced while holding only one record per chunk in memory.

Both premises are asserted rather than trusted. Every VCF stream is checked to
be non-decreasing as it is consumed, and every ``.mtx`` is checked to be
row-major. A tripped assertion means cellsnp-lite did not emit in target-file
order, and this script fails loudly rather than silently writing a mis-ordered
directory.

The invariant that makes the result usable: matrix row *i* must correspond to
record *i* of ``cellSNP.base.vcf`` and record *i* of ``cellSNP.cells.vcf``.
vireo's ``read_cellSNP`` pairs matrix rows with VCF records by position only (it
does not sort), so one permutation is applied to all five files -- the two VCFs
are written in merged order, and the three matrices have their row indices
rewritten through that same permutation. The matrices are emitted row-major in
the merged order, so they stay aligned with the VCFs line for line.

Output filenames match what vireo's ``read_cellSNP`` expects, and the
compression of each output matches the corresponding per-chunk file (whatever
``cellsnp-lite --gzip`` produced).

Usage:
    merge_cellsnp.py <out_dir> <contig_order.txt> <chunk_dir> [<chunk_dir> ...]
"""

import gzip
import heapq
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


def read_contig_ranks(path):
    """Map contig name -> rank from the sorted target VCF's contig order.

    Taken from the order contigs actually appear in the sorted VCF rather than
    from its ``##contig`` header lines, which need not be present, complete, or
    in agreement with the records.
    """
    ranks = {}
    with open(path) as handle:
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


class VcfPairStream:
    """Forward-only reader over one chunk's base.vcf and cells.vcf in lockstep.

    Record *i* of one file pairs with record *i* of the other, so the two are
    always advanced together and never independently. Coordinate order is
    verified as records are consumed.
    """

    def __init__(self, chunk_dir, contig_ranks):
        self.chunk_dir = chunk_dir
        self.contig_ranks = contig_ranks
        self.base_handle = open_maybe_gz(resolve(chunk_dir, BASE_VCF))
        self.cells_handle = open_maybe_gz(resolve(chunk_dir, CELLS_VCF))
        self.base_header, base_first = self._split_header(self.base_handle)
        self.cells_header, cells_first = self._split_header(self.cells_handle)
        self.base_records = self._records(self.base_handle, base_first)
        self.cells_records = self._records(self.cells_handle, cells_first)
        self.record_count = 0
        self.previous_key = None
        self.key = None
        self.base_line = None
        self.cells_line = None
        self.advance()

    @staticmethod
    def _split_header(handle):
        """Consume the leading ``#`` lines, returning them and the first record."""
        header = []
        for line in handle:
            if line.startswith("#"):
                header.append(line)
            elif line.strip():
                return header, line
        return header, None

    @staticmethod
    def _records(handle, first_line):
        if first_line is not None:
            yield first_line
        for line in handle:
            if line.strip():
                yield line

    def advance(self):
        """Load the next paired record, or set key to None at end of stream."""
        base_line = next(self.base_records, None)
        cells_line = next(self.cells_records, None)

        if (base_line is None) != (cells_line is None):
            sys.exit(
                f"{self.chunk_dir}: base and cells VCFs have different record "
                f"counts; they must be row-for-row parallel"
            )
        if base_line is None:
            self.key = None
            self.base_line = None
            self.cells_line = None
            return

        chrom, pos = base_line.split("\t", 2)[:2]
        rank = self.contig_ranks.get(chrom)
        if rank is None:
            sys.exit(
                f"{self.chunk_dir}: contig {chrom!r} is absent from the target "
                f"VCF contig order; this chunk was not produced from that VCF"
            )

        key = (rank, int(pos))
        if self.previous_key is not None and key < self.previous_key:
            sys.exit(
                f"{self.chunk_dir}: records are not in coordinate order "
                f"({self.previous_key} then {key}). cellsnp-lite did not emit "
                f"sites in target-file order, so a streaming merge cannot "
                f"restore the coordinate ordering"
            )

        self.previous_key = key
        self.key = key
        self.base_line = base_line
        self.cells_line = cells_line
        self.record_count += 1

    def close(self):
        self.base_handle.close()
        self.cells_handle.close()


def merge_vcfs(chunk_dirs, contig_ranks, out_dir):
    """k-way merge both VCFs into coordinate order.

    Returns ``(row_map, total_records, base_name, cells_name)`` where
    ``row_map[chunk][local_row - 1]`` is the 1-based row that chunk's record
    occupies in the merged output. That map is what the matrices are rewritten
    through, so the permutation applied here is applied identically there.
    """
    streams = [VcfPairStream(d, contig_ranks) for d in chunk_dirs]
    base_name = basename_of(chunk_dirs[0], BASE_VCF)
    cells_name = basename_of(chunk_dirs[0], CELLS_VCF)
    open_base = gzip.open if base_name.endswith(".gz") else open
    open_cells = gzip.open if cells_name.endswith(".gz") else open

    row_map = [[] for _ in streams]

    # Ties (co-located records split across chunks) break on chunk index, so
    # the merged order is fully determined by the inputs.
    heap = [(s.key, i) for i, s in enumerate(streams) if s.key is not None]
    heapq.heapify(heap)

    total_records = 0
    with open_base(os.path.join(out_dir, base_name), "wt") as base_out, \
            open_cells(os.path.join(out_dir, cells_name), "wt") as cells_out:
        # Headers come from chunk 0; every chunk saw the same BAM and barcodes.
        base_out.writelines(streams[0].base_header)
        cells_out.writelines(streams[0].cells_header)

        while heap:
            _, chunk_index = heapq.heappop(heap)
            stream = streams[chunk_index]

            base_out.write(stream.base_line)
            cells_out.write(stream.cells_line)
            total_records += 1
            row_map[chunk_index].append(total_records)

            stream.advance()
            if stream.key is not None:
                heapq.heappush(heap, (stream.key, chunk_index))

    for stream in streams:
        stream.close()

    return row_map, total_records, base_name, cells_name


class MtxRowStream:
    """Forward-only reader over one chunk's .mtx, yielding entries a row at a time.

    Assumes (and verifies) row-major entry order, which is what lets the merge
    pull rows in ascending local order without ever seeking backwards.
    """

    def __init__(self, path):
        self.path = path
        self.handle = open_maybe_gz(path)
        self.header = []
        self.dims = None
        self.pending = None

        for line in self.handle:
            if line.startswith("%"):
                self.header.append(line)
                continue
            if not line.strip():
                continue
            parts = line.split()
            self.dims = (int(parts[0]), int(parts[1]), int(parts[2]))
            break

        if self.dims is None:
            sys.exit(f"No MatrixMarket dimension line found in {path}")

        self._advance()

    def _advance(self):
        for line in self.handle:
            if not line.strip():
                continue
            parts = line.split()
            self.pending = (int(parts[0]), parts[1], parts[2])
            return
        self.pending = None

    def entries_for_row(self, row):
        """Yield ``(col, value)`` for ``row``; rows must be requested ascending."""
        while self.pending is not None and self.pending[0] <= row:
            if self.pending[0] < row:
                sys.exit(
                    f"{self.path}: entries are not row-major (row "
                    f"{self.pending[0]} appears after row {row}); the merge "
                    f"reads each chunk forward-only and cannot reorder them"
                )
            yield self.pending[1], self.pending[2]
            self._advance()

    def close(self):
        self.handle.close()


def merge_mtx(chunk_dirs, candidates, out_dir, row_map, total_rows):
    """Rewrite one .mtx tag's row indices through the merged-VCF permutation."""
    out_name = basename_of(chunk_dirs[0], candidates)
    out_path = os.path.join(out_dir, out_name)
    open_out = gzip.open if out_name.endswith(".gz") else open

    # Dimension pass: the output header needs the totals before any entry can
    # be written, so collect them before opening the emit streams.
    header_lines = None
    ncols = None
    total_nnz = 0
    for chunk_index, chunk_dir in enumerate(chunk_dirs):
        stream = MtxRowStream(resolve(chunk_dir, candidates))
        nrows, chunk_ncols, nnz = stream.dims
        stream.close()

        if header_lines is None:
            header_lines = stream.header
        if ncols is None:
            ncols = chunk_ncols
        elif chunk_ncols != ncols:
            sys.exit(
                f"{out_name}: column count mismatch ({chunk_ncols} in "
                f"{chunk_dir} vs {ncols} in {chunk_dirs[0]}); barcode sets differ"
            )
        if nrows != len(row_map[chunk_index]):
            sys.exit(
                f"{out_name}: {chunk_dir} declares {nrows} rows but its VCFs "
                f"hold {len(row_map[chunk_index])} records"
            )
        total_nnz += nnz

    # Invert the permutation so the merged rows can be emitted in order. Within
    # any one chunk the local rows are still visited ascending, which is what
    # keeps each source stream forward-only.
    origin = [None] * (total_rows + 1)
    for chunk_index, global_rows in enumerate(row_map):
        for local_row, global_row in enumerate(global_rows, start=1):
            origin[global_row] = (chunk_index, local_row)

    streams = [MtxRowStream(resolve(d, candidates)) for d in chunk_dirs]
    written_nnz = 0
    with open_out(out_path, "wt") as out:
        out.writelines(header_lines)
        out.write(f"{total_rows} {ncols} {total_nnz}\n")

        for global_row in range(1, total_rows + 1):
            chunk_index, local_row = origin[global_row]
            for col, value in streams[chunk_index].entries_for_row(local_row):
                out.write(f"{global_row} {col} {value}\n")
                written_nnz += 1

    for stream in streams:
        stream.close()

    if written_nnz != total_nnz:
        sys.exit(
            f"{out_name}: wrote {written_nnz} entries but the chunk headers "
            f"declared {total_nnz}"
        )

    return out_name, total_rows, ncols, total_nnz


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)

    out_dir = sys.argv[1]
    contig_order_path = sys.argv[2]
    # Sorted only so that ties between co-located records break deterministically;
    # the output order itself now comes from the coordinates, not the chunk order.
    chunk_dirs = sorted(sys.argv[3:], key=os.path.basename)

    os.makedirs(out_dir, exist_ok=True)
    contig_ranks = read_contig_ranks(contig_order_path)

    n_cells = merge_samples(chunk_dirs, out_dir)

    row_map, total_records, base_name, cells_name = merge_vcfs(
        chunk_dirs, contig_ranks, out_dir
    )
    if total_records == 0:
        sys.exit(f"No records found across {len(chunk_dirs)} chunk directories")
    print(f"  {base_name} / {cells_name}: {total_records} SNPs", file=sys.stderr)

    for tag, candidates in MTX_FILES.items():
        mtx_name, mtx_rows, mtx_cols, mtx_nnz = merge_mtx(
            chunk_dirs, candidates, out_dir, row_map, total_records
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
        f"{total_records} SNPs x {n_cells} cells, coordinate-sorted",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
