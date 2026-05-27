#!/usr/bin/env python3
"""
IsoCAPE: Site Annotator (Streaming, GTF-anchored)
--------------------------------------------------
Detects cryptic and alternative polyadenylation sites
that are NOT recorded in GTF annotation.

Decision logic per read:
  1. UNK_GENE                    → skip (no gene context)
  2. known 3'end ±known_window   → known (default 50bp, matches read scatter)
  3. no PAS signal               → skip (noise)
  4. intronic + PAS              → CE (not affected by known_window)
  5. genic + PAS + dist > known_window from any known → PA (truly unannotated)

--known-window 50bp : read scatter window for known 3'end matching.
  Covers ~20bp biological cleavage variation + ~5bp alignment noise.
  CE detection is unaffected — intronic reads are CE regardless of
  distance to known sites.
  PA = genuinely unannotated, > known_window from any GTF 3'end.



Validation (separate workflow):
    IsoCAPE known site:  GENE_known_{coord}
    IsoDecipher G-group: GENE_G0, rep_coord ≈ coord
    → coord matching → Pearson correlation (target r > 0.99)

Usage:
python isocape/annotator/site_annotator.py \
    --parquet path/to/reads.parquet \
    --gtf     path/to/Homo_sapiens.GRCh38.115.gtf \
    --db      path/to/Homo_sapiens.GRCh38.115.gtf.db \
    --ref     path/to/ref.fa \
    --out     path/to/cryptic_sites.parquet \

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
PAS_WINDOW    = 60  # bp upstream of ref_end to scan for PAS signal
                    # covers ~20bp cleavage variation + ~5bp alignment noise
                    # + 30bp typical PAS-to-cleavage distance
PAS_MIN_DIST  = 10  # minimum distance from ref_end to PAS signal
                    # PAS normally 10-30bp upstream of cleavage site
                    # signals closer than 10bp are likely false positives


# ---------------------------------------------------------------------------
# Site name registry
# ---------------------------------------------------------------------------

class SiteRegistry:
    """
    Assigns stable unique names across batches.
    Intronic sites  → GENE_CE1, GENE_CE2 ...
    Genic sites     → GENE_PA1, GENE_PA2 ...
    Known 3' ends   → GENE_known_1, GENE_known_2 ... (--include-known)
    """

    def __init__(self):
        self._ce_counts    = defaultdict(int)
        self._pa_counts    = defaultdict(int)
        self._known_counts = defaultdict(int)
        self._pos_to_name  = {}

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

    def get_known(self, gene, chrom, pos, strand):
        key = (gene, chrom, pos, strand, 'known')
        if key not in self._pos_to_name:
            self._known_counts[gene] += 1
            self._pos_to_name[key] = f"{gene}_known_{self._known_counts[gene]}"
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

    @property
    def n_known(self):
        return sum(self._known_counts.values())


# ---------------------------------------------------------------------------
# PAS signal check
# ---------------------------------------------------------------------------

def check_pas_signal(fasta, chrom, ref_end, strand, window=PAS_WINDOW,
                     min_dist=PAS_MIN_DIST):
    """
    Scan upstream of ref_end for PAS signal.

    Parameters
    ----------
    window   : bp upstream to scan (default 60bp)
    min_dist : minimum distance from ref_end to PAS (default 10bp)
               PAS normally 10-30bp upstream of cleavage site.
               Signals closer than min_dist are likely false positives.
    """
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

        # upstream[-1] is closest to ref_end, upstream[0] is furthest
        # min_dist: PAS must be at least min_dist bp from ref_end
        # i.e. PAS must be in upstream[0 : len(upstream) - min_dist]
        search_region = upstream[:len(upstream) - min_dist]

        for signal in PAS_SIGNALS:
            if signal in search_region:
                return signal
        return None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Annotate one batch
# ---------------------------------------------------------------------------

def annotate_batch(df, parser, fasta, registry, known_window=50):
    """
    Returns DataFrame: known + CE + PA.

    known_window: distance from GTF 3'end for known site matching (default 50bp)
                  covers 90bp read positional scatter
    CE is window-independent — intronic reads are CE regardless of
    distance to known sites.

    PA = genic + PAS + distance > known_window from any known 3'end
    """
    import pandas as pd

    rows = []

    for row in df.itertuples(index=False):
        chrom   = row.chrom
        ref_end = row.ref_end
        strand  = row.strand
        gn      = row.gn

        # 1. Skip UNK_GENE
        if gn == 'UNK_GENE':
            continue

        # 2. Known 3'end check — use known_window (read scatter)
        known_hits = parser.query_3p_end(chrom, ref_end, strand, window=known_window)
        if known_hits:
            gene    = known_hits[0][0]
            site_id = registry.get_known(gene, chrom, ref_end, strand)
            rows.append({
                'cb':        row.cb,
                'ub':        row.ub,
                'gn':        gn,
                'chrom':     chrom,
                'ref_end':   ref_end,
                'strand':    strand,
                'site_id':   site_id,
                'site_type': 'known',
                'pas_signal': None,
                'gene':      gene,
            })
            continue

        # 3. PAS signal required for CE and PA
        pas = check_pas_signal(fasta, chrom, ref_end, strand)
        if not pas:
            continue

        # 4. Classify: intronic (CE) or genic (PA)
        # CE requires:
        #   a) in protein_coding intron (query_intron hit)
        #   b) NOT in any exon (query_exon check)
        #   c) GN tag matches intron gene (consistency check)
        #      → prevents neighboring gene reads being mis-assigned
        intron_hits = parser.query_intron(chrom, ref_end, strand)
        if intron_hits:
            exon_hits = parser.query_exon(chrom, ref_end, strand)
            if not exon_hits:
                intron_gene = intron_hits[0][0]
                # GN tag must match intron gene — else skip (neighboring gene noise)
                if gn == intron_gene:
                    site_id   = registry.get_ce(intron_gene, chrom, ref_end, strand)
                    site_type = 'CE'
                    gene      = intron_gene
                else:
                    continue  # GN tag mismatch → skip
            else:
                # In exon of some transcript → PA
                gene      = gn
                site_id   = registry.get_pa(gene, chrom, ref_end, strand)
                site_type = 'PA'
        else:
            gene      = gn
            site_id   = registry.get_pa(gene, chrom, ref_end, strand)
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
    parser.add_argument("--known-window", type=int, default=50,
                        help="Read scatter window for known 3' end matching (default: 50bp). "
                             "Reads within this distance of a GTF annotated 3'end are labeled "
                             "'known'. Should match read length (~90bp for 10x v3). "
                             "CE detection is unaffected by this parameter.")
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

    print(f"[IsoCAPE] Streaming (batch={args.batch_size:,}, known_window=±{args.known_window}bp)")
    print(f"[IsoCAPE] Output: known + CE + PA")
    print(f"[IsoCAPE] Skipping: UNK_GENE | no PAS signal (for CE/PA)")

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

        df_out = annotate_batch(
            df, gtf_parser, fasta, registry,
            known_window=args.known_window
        )

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
              f"known={registry.n_known:,} "
              f"CE={registry.n_ce:,} PA={registry.n_pa:,} | "
              f"unique={registry.n_unique:,}")

    if writer:
        writer.close()

    fasta.close()

    print(f"\n[SUCCESS] IsoCAPE complete.")
    print(f"  Total reads:        {total_reads:,}")
    print(f"  VALID_PAS reads:    {total_valid:,}")
    print(f"  IsoCAPE site reads: {total_sites:,}")
    print(f"  Unique known sites: {registry.n_known:,}")
    print(f"  Unique CE sites:    {registry.n_ce:,}")
    print(f"  Unique PA sites:    {registry.n_pa:,}")
    print(f"  Total unique sites: {registry.n_unique:,}")
    print(f"  IsoCAPE rate:       {total_sites/total_valid:.3%} of VALID_PAS")
    print(f"  Output:             {args.out}")


if __name__ == "__main__":
    main()
