#!/usr/bin/env python3
"""
IsoCAPE: Site Annotator (Streaming, GTF-anchored)
--------------------------------------------------
Detects cryptic and alternative polyadenylation sites
that are NOT recorded in GTF annotation.

Decision logic per read:
  1. UNK_GENE          → skip (no gene context)
  2. known 3'end ±10bp → skip (IsoDecipher handles these)
  3. no PAS signal     → skip (noise)
  4. intronic + PAS    → GENE_CE (cryptic/premature polyadenylation)
  5. genic + PAS       → GENE_PA (alternative polyadenylation,
                          shorter OR longer than annotated 3'end)

window default = 10bp to match IsoDecipher's tolerance.

Downstream usage:
    apa     = sc.read_h5ad("isodecipher_output.h5ad")   # known APA
    cryptic = sc.read_h5ad("isocape_output.h5ad")        # CE + PA sites
    combined = ad.concat([apa, cryptic], axis=1)

Usage:
python isocape/annotator/site_annotator.py \
    --parquet path/to/reads.parquet \
    --gtf     path/to/Homo_sapiens.GRCh38.115.gtf \
    --db      path/to/Homo_sapiens.GRCh38.115.gtf.db \
    --ref     path/to/ref.fa \
    --out     path/to/cryptic_sites.parquet \
    [--window 10] \
    [--batch-size 500000]
"""

import argparse
import sys
import os
from collections import defaultdict

import pyarrow as pa
import pyarrow.parquet as pq
import pysam

sys.path.insert(0, os.path.dirname(__file__))
from gtf_parser import GTFParser


# ---------------------------------------------------------------------------
# PAS signal variants
# ---------------------------------------------------------------------------
PAS_SIGNALS = {
    "AATAAA", "ATTAAA", "AGTAAA", "TATAAA",
    "CATAAA", "GATAAA", "AATATA", "AATACA",
    "AATAGA", "AATGAA", "ACTAAA", "AACAAA",
}
PAS_WINDOW = 40


# ---------------------------------------------------------------------------
# Site name registry
# ---------------------------------------------------------------------------

class SiteRegistry:
    """
    Assigns stable unique names across batches.
    Intronic sites  → GENE_CE1, GENE_CE2 ...
    Genic sites     → GENE_PA1, GENE_PA2 ...
    """

    def __init__(self):
        self._ce_counts   = defaultdict(int)
        self._pa_counts   = defaultdict(int)
        self._pos_to_name = {}

    def get_ce(self, gene, chrom, pos, strand):
        key = (gene, chrom, pos, strand, 'ce')
        if key not in self._pos_to_name:
            self._ce_counts[gene] += 1
            self._pos_to_name[key] = f"{gene}_CE{self._ce_counts[gene]}"
        return self._pos_to_name[key]

    def get_pa(self, gene, chrom, pos, strand):
        key = (gene, chrom, pos, strand, 'pa')
        if key not in self._pos_to_name:
            self._pa_counts[gene] += 1
            self._pos_to_name[key] = f"{gene}_PA{self._pa_counts[gene]}"
        return self._pos_to_name[key]

    @property
    def n_unique(self):
        return len(self._pos_to_name)

    @property
    def n_ce(self):
        return sum(self._ce_counts.values())

    @property
    def n_pa(self):
        return sum(self._pa_counts.values())


# ---------------------------------------------------------------------------
# PAS signal check
# ---------------------------------------------------------------------------

def check_pas_signal(fasta, chrom, ref_end, strand, window=PAS_WINDOW):
    try:
        chrom_len = fasta.get_reference_length(chrom)
        if strand == '+':
            fetch_start = max(ref_end - window, 0)
            fetch_end   = ref_end
            upstream    = fasta.fetch(chrom, fetch_start, fetch_end).upper()
        else:
            fetch_start = ref_end
            fetch_end   = min(ref_end + window, chrom_len)
            upstream    = fasta.fetch(chrom, fetch_start, fetch_end).upper()
            upstream    = upstream[::-1]

        for signal in PAS_SIGNALS:
            if signal in upstream:
                return signal
        return None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Annotate one batch
# ---------------------------------------------------------------------------

def annotate_batch(df, parser, fasta, registry, window):
    """
    Returns DataFrame of IsoCAPE sites (CE + PA) only.
    Skips: UNK_GENE, known 3'ends, no-PAS reads.
    """
    import pandas as pd

    rows = []

    for row in df.itertuples(index=False):
        chrom   = row.chrom
        ref_end = row.ref_end
        strand  = row.strand
        gn      = row.gn

        # 1. Skip UNK_GENE — no gene context
        if gn == 'UNK_GENE':
            continue

        # 2. Skip known 3'ends (±window = IsoDecipher's territory)
        if parser.query_3p_end(chrom, ref_end, strand, window=window):
            continue

        # 3. PAS signal required
        pas = check_pas_signal(fasta, chrom, ref_end, strand)
        if not pas:
            continue

        # 4. Classify: intronic (CE) or genic (PA)
        intron_hits = parser.query_intron(chrom, ref_end, strand)
        if intron_hits:
            gene    = intron_hits[0][0]
            site_id = registry.get_ce(gene, chrom, ref_end, strand)
            site_type = 'CE'
        else:
            gene    = gn  # use Cell Ranger GN tag
            site_id = registry.get_pa(gene, chrom, ref_end, strand)
            site_type = 'PA'

        rows.append({
            'cb':        row.cb,
            'ub':        row.ub,
            'gn':        gn,
            'chrom':     chrom,
            'ref_end':   ref_end,
            'strand':    strand,
            'site_id':   site_id,
            'site_type': site_type,
            'pas_signal': pas,
            'gene':      gene,
        })

    if not rows:
        return None
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Output schema
# ---------------------------------------------------------------------------

OUT_SCHEMA = pa.schema([
    pa.field("cb",         pa.string()),
    pa.field("ub",         pa.string()),
    pa.field("gn",         pa.string()),
    pa.field("chrom",      pa.string()),
    pa.field("ref_end",    pa.int64()),
    pa.field("strand",     pa.string()),
    pa.field("site_id",    pa.string()),
    pa.field("site_type",  pa.string()),
    pa.field("pas_signal", pa.string()),
    pa.field("gene",       pa.string()),
])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="IsoCAPE Site Annotator — cryptic & alternative polyA"
    )
    parser.add_argument("--parquet",    required=True)
    parser.add_argument("--gtf",        required=True)
    parser.add_argument("--db",         default=None)
    parser.add_argument("--ref",        required=True)
    parser.add_argument("--out",        required=True)
    parser.add_argument("--window",     type=int, default=10,
                        help="bp tolerance for known 3' end match "
                             "(default: 10, matches IsoDecipher tolerance)")
    parser.add_argument("--batch-size", type=int, default=500_000)
    parser.add_argument("--include-no-ref", action="store_true",
                        help="Include NO_REF reads (chrom not in reference FASTA). "
                             "Default: only VALID_PAS reads. "
                             "Use during testing with partial reference.")
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"[IsoCAPE] Building GTF index ...")
    gtf_parser = GTFParser(gtf_path=args.gtf, db_path=args.db)
    gtf_parser.build()

    fasta    = pysam.FastaFile(args.ref)
    registry = SiteRegistry()
    pf       = pq.ParquetFile(args.parquet)
    writer   = None

    # Determine which priming labels to include
    if args.include_no_ref:
        valid_labels = {'VALID_PAS', 'NO_REF'}
        print(f"[IsoCAPE] Including NO_REF reads (testing mode with partial reference)")
    else:
        valid_labels = {'VALID_PAS'}
        print(f"[IsoCAPE] VALID_PAS only (use --include-no-ref for partial reference)")

    total_reads = 0
    total_valid = 0
    total_sites = 0
    batch_num   = 0

    print(f"[IsoCAPE] Streaming (batch={args.batch_size:,}, window=±{args.window}bp)")
    print(f"[IsoCAPE] Skipping: UNK_GENE | known 3'ends | no PAS signal")
    print(f"[IsoCAPE] Keeping:  intronic+PAS → CE | genic+PAS → PA")

    for batch in pf.iter_batches(batch_size=args.batch_size):
        import pandas as pd
        df = batch.to_pandas()
        batch_num   += 1
        total_reads += len(df)

        df = df[df['priming_label'].isin(valid_labels)].copy()
        total_valid += len(df)

        if df.empty:
            continue

        df = df.drop(columns=['sequence'], errors='ignore')

        df_out = annotate_batch(df, gtf_parser, fasta, registry, args.window)

        if df_out is None or df_out.empty:
            continue

        total_sites += len(df_out)

        table = pa.Table.from_pandas(df_out, schema=OUT_SCHEMA, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(args.out, OUT_SCHEMA)
        writer.write_table(table)

        print(f"  [Batch {batch_num}] "
              f"reads={total_reads:,} | "
              f"valid={total_valid:,} | "
              f"sites={total_sites:,} | "
              f"CE={registry.n_ce:,} PA={registry.n_pa:,} | "
              f"unique={registry.n_unique:,}")

    if writer:
        writer.close()

    fasta.close()

    print(f"\n[SUCCESS] IsoCAPE complete.")
    print(f"  Total reads:        {total_reads:,}")
    print(f"  VALID_PAS reads:    {total_valid:,}")
    print(f"  IsoCAPE site reads: {total_sites:,}")
    print(f"  Unique CE sites:    {registry.n_ce:,}")
    print(f"  Unique PA sites:    {registry.n_pa:,}")
    print(f"  Total unique sites: {registry.n_unique:,}")
    print(f"  IsoCAPE rate:       {total_sites/total_valid:.3%} of VALID_PAS")
    print(f"  Output:             {args.out}")


if __name__ == "__main__":
    main()
