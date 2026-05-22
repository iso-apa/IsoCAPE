#!/usr/bin/env python3
"""
IsoCAPE: BAM to Parquet ETL (Parallel Streaming version)
----------------------------------------------------------
Parallelized by chromosome. Each worker streams reads directly
to disk using PyArrow's incremental ParquetWriter — RAM usage
is O(1) per worker regardless of chromosome size.

Merge step also uses PyArrow streaming — reads one chromosome
parquet at a time, writes to final output, deletes from RAM.

Usage:
python isocape/scripts/bam_to_parquet_parallel.py \
    --bam path/to/input.bam \
    --ref path/to/ref.fa \
    --out path/to/output.parquet \
    --barcodes path/to/filtered_barcodes.csv \
    [--cores 4] \
    [--window 80] \
    [--min-phred 20] \
    [--priming-window 20] \
    [--priming-threshold 8] \
    [--chroms chr1,chr2,chr7] \
    [--chunk-size 10000]
"""

import argparse
import gzip
import os
import tempfile
from multiprocessing import Pool, cpu_count

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pandas as pd
import pysam
from collections import defaultdict


# Standard chromosomes — skip contigs by default
STANDARD_CHROMS = set(
    [f"chr{i}" for i in range(1, 23)] +
    ["chrX", "chrY", "chrM"] +
    [str(i) for i in range(1, 23)] +
    ["X", "Y", "MT"]
)

# PyArrow schema
SCHEMA = pa.schema([
    pa.field("cb",            pa.string()),
    pa.field("ub",            pa.string()),
    pa.field("gn",            pa.string()),
    pa.field("sequence",      pa.string()),
    pa.field("strand",        pa.string()),
    pa.field("chrom",         pa.string()),
    pa.field("ref_end",       pa.int64()),
    pa.field("priming_label", pa.string()),
])


# ---------------------------------------------------------------------------
# Barcode loader
# ---------------------------------------------------------------------------

def detect_sep(path):
    p = path.lower().replace('.gz', '')
    return '\t' if p.endswith('.tsv') else ','


def load_barcodes(barcodes_path, barcode_col=None, sep=None):
    if sep is None:
        sep = detect_sep(barcodes_path)

    opener = gzip.open if barcodes_path.endswith('.gz') else open
    with opener(barcodes_path, 'rt') as fh:
        first_line = fh.readline().strip()

    fields = first_line.split(sep)
    if barcode_col is None:
        barcode_col = 0
        for i, f in enumerate(fields):
            cleaned = f.replace('-1', '').split('-')[0]
            if all(c in 'ACGTacgt' for c in cleaned) and len(cleaned) >= 12:
                barcode_col = i
                break

    df = pd.read_csv(barcodes_path, header=None, sep=sep)
    barcodes = set(df[barcode_col].astype(str).str.replace(r'-\d+$', '', regex=True))
    print(f"[ETL] Loaded {len(barcodes):,} valid cell barcodes")
    return barcodes


# ---------------------------------------------------------------------------
# Internal priming check
# ---------------------------------------------------------------------------

def check_internal_priming(fasta, chrom, ref_end, strand, window, threshold):
    """
    Check downstream sequence for A-tract (internal priming signal).

    Returns:
        'VALID_PAS'      — ref exists, no A-tract downstream
        'INTERNAL_PRIME' — ref exists, A-tract found → likely internal priming
        'NO_REF'         — chrom not in reference FASTA → cannot verify
                           (use partial reference during testing;
                            re-run with full hg38 for production)
    """
    try:
        # Check chrom exists in reference before fetching
        if chrom not in fasta.references:
            return 'NO_REF'

        chrom_len = fasta.get_reference_length(chrom)
        if strand == '+':
            fetch_start = ref_end
            fetch_end   = min(ref_end + window, chrom_len)
            downstream  = fasta.fetch(chrom, fetch_start, fetch_end).upper()
        else:
            fetch_start = max(ref_end - window, 0)
            fetch_end   = ref_end
            downstream  = fasta.fetch(chrom, fetch_start, fetch_end).upper()
            downstream  = downstream[::-1]

        max_run = current_run = 0
        for base in downstream:
            if base == 'A':
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 0

        return 'INTERNAL_PRIME' if max_run >= threshold else 'VALID_PAS'
    except Exception:
        return 'NO_REF'


# ---------------------------------------------------------------------------
# Per-chromosome worker — streaming
# ---------------------------------------------------------------------------

def process_chromosome(args):
    (
        bam_path, ref_path, chrom,
        valid_barcodes, min_phred,
        priming_window, priming_threshold,
        window_size, chunk_size, tmp_dir
    ) = args

    bam   = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(ref_path)

    tmp_path = os.path.join(tmp_dir, f"{chrom}.parquet")
    writer   = None

    seen_umis  = set()
    buffer     = defaultdict(list)
    n_retained = 0

    def flush_buffer():
        nonlocal writer
        if not buffer["cb"]:
            return
        table = pa.table(buffer, schema=SCHEMA)
        if writer is None:
            writer = pq.ParquetWriter(tmp_path, SCHEMA)
        writer.write_table(table)
        for key in buffer:
            buffer[key].clear()

    for read in bam.fetch(chrom):
        if read.is_unmapped:
            continue
        if not (read.has_tag("CB") and read.has_tag("UB")):
            continue

        cb = read.get_tag("CB")
        ub = read.get_tag("UB")

        if valid_barcodes is not None:
            cb_stripped = cb.replace('-1', '').split('-')[0]
            if cb_stripped not in valid_barcodes:
                continue

        seq = read.query_sequence
        if seq is None or len(seq) < window_size:
            continue
        window_seq = seq[-window_size:]

        quals = read.query_qualities
        if quals is not None:
            if np.mean(quals[-window_size:]) < min_phred:
                continue

        umi_key = (cb, ub)
        if umi_key in seen_umis:
            continue
        seen_umis.add(umi_key)

        strand  = '-' if read.is_reverse else '+'
        ref_end = read.reference_end

        priming_label = check_internal_priming(
            fasta, chrom, ref_end, strand,
            priming_window, priming_threshold
        )

        buffer["cb"].append(cb)
        buffer["ub"].append(ub)
        buffer["gn"].append(read.get_tag("GN") if read.has_tag("GN") else "UNK_GENE")
        buffer["sequence"].append(window_seq)
        buffer["strand"].append(strand)
        buffer["chrom"].append(chrom)
        buffer["ref_end"].append(ref_end)
        buffer["priming_label"].append(priming_label)
        n_retained += 1

        if n_retained % chunk_size == 0:
            flush_buffer()

    flush_buffer()

    if writer:
        writer.close()

    bam.close()
    fasta.close()

    if n_retained > 0:
        print(f"  [Worker] {chrom}: {n_retained:,} reads retained")
        return tmp_path
    else:
        print(f"  [Worker] {chrom}: 0 reads (skipped)")
        return None


# ---------------------------------------------------------------------------
# Streaming merge — PyArrow only, no pandas concat
# ---------------------------------------------------------------------------

def streaming_merge(tmp_paths, output_path):
    """
    Merge chromosome parquets into one output file.
    Reads one file at a time → writes → deletes from RAM.
    RAM usage: O(one chromosome) not O(all chromosomes).
    """
    print(f"[ETL] Streaming merge → {output_path}")

    writer = None
    total  = 0

    for path in tmp_paths:
        if path is None or not os.path.exists(path):
            continue

        table = pq.read_table(path)          # read one chrom
        n = len(table)

        if writer is None:
            writer = pq.ParquetWriter(output_path, SCHEMA)

        writer.write_table(table)            # write to output
        total += n
        del table                            # free RAM immediately

    if writer:
        writer.close()

    print(f"[ETL] Merge complete: {total:,} total reads")
    return total


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="IsoCAPE BAM to Parquet ETL (Parallel Streaming)"
    )
    parser.add_argument("--bam",               required=True)
    parser.add_argument("--ref",               required=True)
    parser.add_argument("--out",               required=True)
    parser.add_argument("--barcodes",          default=None)
    parser.add_argument("--barcode-col",       type=int, default=None)
    parser.add_argument("--barcode-sep",       default=None)
    parser.add_argument("--cores",             type=int, default=None)
    parser.add_argument("--chroms",            default=None,
                        help="Comma-separated chroms. Default: standard chr1-22,X,Y,M. "
                             "Use 'all' for all contigs.")
    parser.add_argument("--chunk-size",        type=int, default=10_000)
    parser.add_argument("--min-phred",         type=int, default=20)
    parser.add_argument("--priming-window",    type=int, default=20)
    parser.add_argument("--priming-threshold", type=int, default=8)
    parser.add_argument("--window",            type=int, default=100)
    return parser.parse_args()


def main():
    args = parse_args()

    valid_barcodes = None
    if args.barcodes:
        valid_barcodes = load_barcodes(
            args.barcodes,
            barcode_col=args.barcode_col,
            sep=args.barcode_sep,
        )

    bam = pysam.AlignmentFile(args.bam, "rb")
    all_chroms = [sq['SN'] for sq in bam.header['SQ']]
    bam.close()
    print(f"[ETL] Found {len(all_chroms)} chromosomes in BAM")

    if args.chroms == 'all':
        chroms = all_chroms
        print(f"[ETL] Using all {len(chroms)} chromosomes")
    elif args.chroms:
        chroms = [c.strip() for c in args.chroms.split(',')]
        print(f"[ETL] Using {len(chroms)} user-specified chromosomes")
    else:
        chroms = [c for c in all_chroms if c in STANDARD_CHROMS]
        print(f"[ETL] Using {len(chroms)} standard chromosomes "
              f"(use --chroms all to include contigs)")

    n_cores = args.cores or cpu_count()
    n_cores = min(n_cores, len(chroms))
    print(f"[ETL] Using {n_cores} cores | chunk_size={args.chunk_size:,}")

    with tempfile.TemporaryDirectory() as tmp_dir:
        worker_args = [
            (
                args.bam, args.ref, chrom,
                valid_barcodes, args.min_phred,
                args.priming_window, args.priming_threshold,
                args.window, args.chunk_size, tmp_dir
            )
            for chrom in chroms
        ]

        print(f"[ETL] Starting parallel streaming extraction ...")
        with Pool(processes=n_cores) as pool:
            tmp_paths = pool.map(process_chromosome, worker_args)

        # Streaming merge — no pandas, no RAM explosion
        total = streaming_merge(tmp_paths, args.out)

    if total == 0:
        print("[ERROR] No reads retained.")
        return

    # Quick summary from output file
    print(f"\n[SUCCESS] Saved {total:,} reads → {args.out}")

    # Read just the metadata for summary (no full load)
    pf = pq.read_table(args.out, columns=["cb", "gn", "priming_label"])
    print(f"  Unique cells (CB): {pf['cb'].n_unique():,}")
    print(f"  Unique genes (GN): {pf['gn'].n_unique():,}")
    valid  = (pf['priming_label'] == 'VALID_PAS').sum().as_py()
    prime  = (pf['priming_label'] == 'INTERNAL_PRIME').sum().as_py()
    no_ref = (pf['priming_label'] == 'NO_REF').sum().as_py()
    unk    = (pf['gn'] == 'UNK_GENE').sum().as_py()
    print(f"  VALID_PAS:         {valid:,}")
    print(f"  INTERNAL_PRIME:    {prime:,}")
    print(f"  NO_REF:            {no_ref:,}  (chrom not in ref — re-run with full hg38)")
    print(f"  UNK_GENE reads:    {unk:,}")
    if no_ref > 0:
        print(f"\n  [NOTE] {no_ref/total:.1%} of reads have NO_REF priming label.")
        print(f"         Internal priming could not be verified for these reads.")
        print(f"         For production use, provide a complete hg38 reference.")
        print(f"         site_annotator.py treats NO_REF reads as unverified —")
        print(f"         use --include-no-ref flag to include them in analysis.")


if __name__ == "__main__":
    main()
