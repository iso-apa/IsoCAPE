#!/usr/bin/env python3
"""
IsoCAPE: BAM to Parquet ETL
-----------------------------
Extracts the final 100bp from each 3' scRNA-seq read,
applies valid barcode filtering, UMI deduplication, Phred QC,
and internal priming labeling. Outputs an AI-ready Parquet file.

No GTF, no gene list, no panel required.

Usage:
python isocape/scripts/bam_to_parquet.py \
    --bam path/to/input.bam \
    --ref path/to/melanoma_ref.fa \
    --out path/to/output.parquet \
    --barcodes path/to/filtered_barcodes.csv \
    [--barcode-col 1] \
    [--barcode-sep ,] \
    [--min-phred 20] \
    [--priming-window 20] \
    [--priming-threshold 8] \
    [--window 100]
"""

import argparse
import gzip
import pysam
import pandas as pd
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="IsoCAPE BAM to Parquet ETL"
    )
    parser.add_argument("--bam",               required=True,  help="Input BAM file")
    parser.add_argument("--ref",               required=True,  help="Reference genome FASTA")
    parser.add_argument("--out",               required=True,  help="Output Parquet file")
    parser.add_argument("--barcodes",          default=None,   help="Valid cell barcodes file (CSV, TSV, or .gz)")
    parser.add_argument("--barcode-col",       type=int, default=None,
                        help="0-based column index containing barcodes. "
                             "Auto-detected if not set.")
    parser.add_argument("--barcode-sep",       default=None,
                        help="Delimiter for barcodes file. Auto-detected from extension if not set.")
    parser.add_argument("--min-phred",         type=int, default=20,  help="Minimum mean Phred score (default: 20)")
    parser.add_argument("--priming-window",    type=int, default=20,  help="Downstream bp to check for A-tract (default: 20)")
    parser.add_argument("--priming-threshold", type=int, default=8,   help="Consecutive A's to flag internal priming (default: 8)")
    parser.add_argument("--window",            type=int, default=100, help="3' end window size in bp (default: 100)")
    return parser.parse_args()


def detect_sep(path):
    """Infer delimiter from file extension."""
    p = path.lower().replace('.gz', '')
    if p.endswith('.tsv'):
        return '\t'
    return ','   # default CSV


def load_barcodes(barcodes_path, barcode_col=None, sep=None):
    """
    Load valid cell barcodes from Cell Ranger output.

    Handles:
    - Single-column TSV/CSV (standard Cell Ranger barcodes.tsv.gz)
    - Two-column CSV with genome prefix (e.g. GRCh38,BARCODE-1)
    - Gzipped or plain text
    - Auto-detects barcode column if not specified

    Parameters
    ----------
    barcodes_path : str
    barcode_col   : int or None — 0-based column index. Auto-detected if None.
    sep           : str or None — delimiter. Auto-detected if None.

    Returns
    -------
    set of barcode strings (suffix -1, -2 etc stripped)
    """
    if sep is None:
        sep = detect_sep(barcodes_path)

    opener = gzip.open if barcodes_path.endswith('.gz') else open

    # Peek at first line to figure out structure
    with opener(barcodes_path, 'rt') as fh:
        first_line = fh.readline().strip()

    fields = first_line.split(sep)
    n_cols = len(fields)

    if barcode_col is None:
        # Heuristic: barcodes look like ACGT...(-N)
        # Find the column that matches a barcode pattern
        barcode_col = 0
        for i, f in enumerate(fields):
            cleaned = f.replace('-1', '').split('-')[0]
            if all(c in 'ACGTacgt' for c in cleaned) and len(cleaned) >= 12:
                barcode_col = i
                break

    print(f"[ETL] Barcode file: {n_cols} column(s), sep='{sep}', "
          f"using column {barcode_col}")

    df = pd.read_csv(barcodes_path, header=None, sep=sep)
    raw = df[barcode_col].astype(str)

    # Strip -1, -2 suffixes (Cell Ranger gem group suffix)
    barcodes = set(raw.str.replace(r'-\d+$', '', regex=True))

    print(f"[ETL] Loaded {len(barcodes):,} valid cell barcodes from {barcodes_path}")
    return barcodes


def check_internal_priming(fasta, chrom, ref_end, strand, window, threshold):
    """
    Query genomic sequence downstream of read 3' end.
    Returns: 'INTERNAL_PRIME' or 'VALID_PAS'
    """
    try:
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

        max_run = 0
        current_run = 0
        for base in downstream:
            if base == 'A':
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 0

        return 'INTERNAL_PRIME' if max_run >= threshold else 'VALID_PAS'

    except Exception as e:
        print(f"[WARN] fasta.fetch failed: {chrom}:{ref_end} strand={strand} ({e})")
        return 'VALID_PAS'


def extract_reads(bam_path, ref_path, min_phred, priming_window,
                  priming_threshold, window_size, valid_barcodes=None):
    """
    Main ETL loop:
    1. Filter reads (CB + UB required, mapped, long enough)
    2. Valid barcode filter (if --barcodes provided)
    3. Phred QC on 3' window
    4. UMI deduplication
    5. Internal priming labeling
    """
    bam   = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(ref_path)

    seen_umis = set()
    rows = []
    stats = defaultdict(int)

    print(f"[ETL] Starting BAM extraction (window={window_size}bp, min_phred={min_phred})")
    if valid_barcodes:
        print(f"[ETL] Barcode filter: ON ({len(valid_barcodes):,} valid barcodes)")
    else:
        print(f"[ETL] Barcode filter: OFF")

    for read in bam.fetch():
        stats['total'] += 1

        if read.is_unmapped:
            stats['unmapped'] += 1
            continue
        if not (read.has_tag("CB") and read.has_tag("UB")):
            stats['no_cb_ub'] += 1
            continue

        cb = read.get_tag("CB")
        ub = read.get_tag("UB")

        # --- Valid barcode filter ---
        if valid_barcodes is not None:
            cb_stripped = cb.replace('-1', '').split('-')[0]
            if cb_stripped not in valid_barcodes:
                stats['invalid_cb'] += 1
                continue

        # --- Sequence extraction ---
        seq = read.query_sequence
        if seq is None or len(seq) < window_size:
            stats['too_short'] += 1
            continue

        window_seq = seq[-window_size:]

        # --- Phred QC on window ---
        quals = read.query_qualities
        if quals is not None:
            window_quals = quals[-window_size:]
            if window_quals.mean() < min_phred:
                stats['low_qual'] += 1
                continue

        # --- UMI deduplication ---
        umi_key = (cb, ub)
        if umi_key in seen_umis:
            stats['umi_dup'] += 1
            continue
        seen_umis.add(umi_key)

        # --- Internal priming label ---
        strand  = '-' if read.is_reverse else '+'
        ref_end = read.reference_end
        chrom   = read.reference_name

        priming_label = check_internal_priming(
            fasta, chrom, ref_end, strand,
            priming_window, priming_threshold
        )
        stats[priming_label] += 1

        rows.append({
            'cb':            cb,
            'ub':            ub,
            'gn':            read.get_tag("GN") if read.has_tag("GN") else "UNK_GENE",
            'sequence':      window_seq,
            'strand':        strand,
            'chrom':         chrom,
            'ref_end':       ref_end,
            'priming_label': priming_label,
        })

        if len(rows) % 500_000 == 0:
            print(f"  [ETL] {len(rows):,} reads retained so far...")

    bam.close()
    fasta.close()

    print(f"\n[ETL] Done.")
    print(f"  Total reads:          {stats['total']:>10,}")
    print(f"  Unmapped:             {stats['unmapped']:>10,}")
    print(f"  No CB/UB:             {stats['no_cb_ub']:>10,}")
    print(f"  Invalid barcode:      {stats['invalid_cb']:>10,}")
    print(f"  Too short (<{window_size}bp):  {stats['too_short']:>10,}")
    print(f"  Low Phred (<{min_phred}):     {stats['low_qual']:>10,}")
    print(f"  UMI duplicates:       {stats['umi_dup']:>10,}")
    print(f"  VALID_PAS retained:   {stats['VALID_PAS']:>10,}")
    print(f"  INTERNAL_PRIME:       {stats['INTERNAL_PRIME']:>10,}")
    print(f"  Total retained:       {len(rows):>10,}")

    return pd.DataFrame(rows)


def main():
    args = parse_args()

    valid_barcodes = None
    if args.barcodes:
        valid_barcodes = load_barcodes(
            args.barcodes,
            barcode_col = args.barcode_col,
            sep         = args.barcode_sep,
        )

    df = extract_reads(
        bam_path          = args.bam,
        ref_path          = args.ref,
        min_phred         = args.min_phred,
        priming_window    = args.priming_window,
        priming_threshold = args.priming_threshold,
        window_size       = args.window,
        valid_barcodes    = valid_barcodes,
    )

    if df.empty:
        print("[ERROR] No reads retained. Check BAM tags and filters.")
        return

    df.to_parquet(args.out, index=False)
    print(f"\n[SUCCESS] Saved {len(df):,} reads → {args.out}")
    print(f"  Unique cells (CB): {df['cb'].nunique():,}")
    print(f"  Unique genes (GN): {df['gn'].nunique():,}")
    print(f"  UNK_GENE reads:    {(df['gn'] == 'UNK_GENE').sum():,}")


if __name__ == "__main__":
    main()
