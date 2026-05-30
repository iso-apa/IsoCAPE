#!/usr/bin/env python3
"""
IsoCAPE: Site Annotator (Streaming, GTF-anchored)
--------------------------------------------------
Detects cryptic and alternative polyadenylation sites
that are NOT recorded in GTF annotation.

Decision logic per read:
  1. UNK_GENE (no PAS, not Alu boundary) → skip
  1b. UNK_GENE + Alu boundary            → AluCE (deep intronic reads)
  2. known 3'end ±known_window           → known (default 50bp)
  3. no PAS signal + not Alu boundary    → skip (noise)
  3b. no PAS signal + Alu boundary       → AluCE
  4. intronic + PAS                      → CE
  5. genic + PAS + dist > known_window   → PA

Alu CE detection:
  Uses RepeatMasker BED file (--alu-bed) when provided.
  Falls back to soft-masking check (hg38 lowercase) otherwise.
  Genuine Alu CE: ref_end is just past Alu element boundary
    → upstream overlaps with Alu element
    → ref_end itself is NOT inside any Alu element
  Internal priming: ref_end is inside Alu body → rejected

Usage:
python isocape/annotator/site_annotator.py \
    --parquet path/to/reads.parquet \
    --gtf     path/to/Homo_sapiens.GRCh38.115.gtf \
    --db      path/to/Homo_sapiens.GRCh38.115.gtf.db \
    --ref     path/to/ref.fa \
    --out     path/to/cryptic_sites.parquet \
    --alu-bed path/to/alu_hg38.bed \
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
PAS_WINDOW    = 60
PAS_MIN_DIST  = 10
ALU_WINDOW    = 50  # bp upstream to scan for Alu overlap


# ---------------------------------------------------------------------------
# Site name registry
# ---------------------------------------------------------------------------

class SiteRegistry:
    def __init__(self):
        self._ce_counts    = defaultdict(int)
        self._pa_counts    = defaultdict(int)
        self._known_counts = defaultdict(int)
        self._alu_counts   = defaultdict(int)
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

    def get_alu(self, gene, chrom, pos, strand):
        key = (gene, chrom, pos, strand, 'alu')
        if key not in self._pos_to_name:
            self._alu_counts[gene] += 1
            self._pos_to_name[key] = f"{gene}_AluCE_{pos}"
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

    @property
    def n_alu(self):
        return sum(self._alu_counts.values())


# ---------------------------------------------------------------------------
# PAS signal check
# ---------------------------------------------------------------------------

def check_pas_signal(fasta, chrom, ref_end, strand, window=PAS_WINDOW,
                     min_dist=PAS_MIN_DIST):
    """Scan upstream of ref_end for PAS signal."""
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

        search_region = upstream[:len(upstream) - min_dist]
        for signal in PAS_SIGNALS:
            if signal in search_region:
                return signal
        return None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Alu element boundary check
# ---------------------------------------------------------------------------

def load_alu_trees(alu_bed_path):
    """
    Load RepeatMasker Alu BED file into per-chromosome interval trees.
    BED format: chrom start end name score strand
    Generated by: gunzip -c rmsk.txt.gz | awk '$13=="Alu"' | ...
    """
    from intervaltree import IntervalTree
    import pandas as pd

    print(f"[IsoCAPE] Loading Alu BED: {alu_bed_path}")
    df = pd.read_csv(alu_bed_path, sep='\t', header=None,
                     names=['chrom','start','end','name','score','strand'])
    trees = defaultdict(IntervalTree)
    for row in df.itertuples(index=False):
        trees[row.chrom][row.start:row.end] = row.name
    print(f"[IsoCAPE] Loaded {len(df):,} Alu elements across {len(trees)} chroms")
    return trees


def is_alu_boundary_rmsk(alu_trees, chrom, ref_end, strand,
                           window=ALU_WINDOW):
    """
    Check Alu CE boundary using RepeatMasker interval tree.

    Genuine Alu CE:
      - upstream window overlaps with an Alu element
      - ref_end itself is NOT inside any Alu element
      → reads terminate just past the Alu 3' end

    Internal priming:
      - ref_end is inside an Alu element
      → reads terminate within the Alu body (A-rich linker)

    Returns
    -------
    (bool, str) : (is_alu_boundary, alu_name or '')
    """
    chrom_key = chrom if chrom.startswith('chr') else f'chr{chrom}'

    # ref_end inside Alu → internal priming, not CE
    in_alu_hits = alu_trees[chrom_key][ref_end - 1:ref_end + 1]
    if in_alu_hits:
        return False, ''

    # upstream overlaps with Alu → at boundary
    if strand == '+':
        up_hits = alu_trees[chrom_key][max(0, ref_end - window):ref_end]
    else:
        up_hits = alu_trees[chrom_key][ref_end:ref_end + window]

    if up_hits:
        alu_name = next(iter(up_hits)).data
        return True, alu_name
    return False, ''


def is_alu_boundary_softmask(fasta, chrom, ref_end, strand,
                               window=ALU_WINDOW, min_alu_pct=0.7):
    """
    Fallback: check Alu CE boundary using hg38 soft-masking (lowercase).
    Less precise than RepeatMasker but requires no extra file.

    Returns
    -------
    (bool, float) : (is_alu_boundary, alu_pct_upstream)
    """
    try:
        chrom_fa  = chrom if chrom.startswith('chr') else f'chr{chrom}'
        chrom_len = fasta.get_reference_length(chrom_fa)
        if strand == '+':
            up_seq = fasta.fetch(chrom_fa, max(0, ref_end - window), ref_end)
            ds_seq = fasta.fetch(chrom_fa, ref_end, min(ref_end + 5, chrom_len))
        else:
            up_seq = fasta.fetch(chrom_fa, ref_end, min(ref_end + window, chrom_len))
            up_seq = up_seq[::-1]
            ds_seq = fasta.fetch(chrom_fa, max(0, ref_end - 5), ref_end)
            ds_seq = ds_seq[::-1]

        if not up_seq:
            return False, 0.0

        n_lower  = sum(1 for c in up_seq if c.islower())
        alu_pct  = n_lower / len(up_seq)
        at_boundary = (alu_pct >= min_alu_pct and
                       len(ds_seq) > 0 and
                       ds_seq == ds_seq.upper() and
                       ds_seq.isalpha())
        return at_boundary, alu_pct
    except Exception:
        return False, 0.0


def check_alu_boundary(fasta, chrom, ref_end, strand, alu_trees=None):
    """
    Unified Alu boundary check.
    Uses RepeatMasker trees if available, falls back to soft-masking.

    Returns
    -------
    (bool, str) : (is_alu_boundary, alu_label)
                  alu_label = Alu family name (rmsk) or 'Alu:{pct}' (softmask)
    """
    if alu_trees is not None:
        is_alu, alu_name = is_alu_boundary_rmsk(alu_trees, chrom, ref_end, strand)
        return is_alu, alu_name
    else:
        is_alu, alu_pct = is_alu_boundary_softmask(fasta, chrom, ref_end, strand)
        return is_alu, f'Alu:{alu_pct:.2f}'


# ---------------------------------------------------------------------------
# Annotate one batch
# ---------------------------------------------------------------------------

def annotate_batch(df, parser, fasta, registry, known_window=50,
                   alu_trees=None):
    """
    Returns DataFrame: known + CE + AluCE + PA.

    Alu CE detection uses RepeatMasker trees (preferred) or soft-masking.
    UNK_GENE reads are allowed through only if at confirmed Alu boundary.
    """
    import pandas as pd

    rows = []

    for row in df.itertuples(index=False):
        chrom   = row.chrom
        ref_end = row.ref_end
        strand  = row.strand
        gn      = row.gn

        # ── 1. UNK_GENE ───────────────────────────────────────────────────
        # Cell Ranger does not assign GN to reads far from exon boundaries.
        # Alu CE reads in deep introns (e.g. MGA intron 17) are often UNK_GENE.
        # Allow through ONLY if at confirmed Alu boundary in a protein_coding intron.
        if gn == 'UNK_GENE':
            # Try read strand first, then opposite strand.
            # Antisense Alu: +strand genes have Alu inserted antisense.
            # oligo-dT captures the Alu A-rich tail on the antisense strand
            # → reads map to - strand (flag=16) even though gene is + strand.
            opposite_strand = '+' if strand == '-' else '-'
            for query_strand in [strand, opposite_strand]:
                intron_hits_unk = parser.query_intron(chrom, ref_end, query_strand)
                if not intron_hits_unk:
                    continue
                exon_hits_unk = parser.query_exon(chrom, ref_end, query_strand)
                if exon_hits_unk:
                    continue
                is_alu, alu_label = check_alu_boundary(
                    fasta, chrom, ref_end, query_strand, alu_trees)
                if is_alu:
                    intron_gene_unk = intron_hits_unk[0][0]
                    site_id = registry.get_alu(
                        intron_gene_unk, chrom, ref_end, strand)
                    rows.append({
                        'cb':        row.cb,
                        'ub':        row.ub,
                        'gn':        gn,
                        'chrom':     chrom,
                        'ref_end':   ref_end,
                        'strand':    strand,
                        'site_id':   site_id,
                        'site_type': 'AluCE',
                        'pas_signal': alu_label,
                        'gene':      intron_gene_unk,
                    })
                    break  # found on this strand, stop
            continue  # UNK_GENE: skip if not Alu boundary

        # ── 2. Known 3'end ────────────────────────────────────────────────
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

        # ── 3. PAS signal check ───────────────────────────────────────────
        pas = check_pas_signal(fasta, chrom, ref_end, strand)

        # ── 3b. No PAS → check for Alu CE ────────────────────────────────
        # Alu CE lacks canonical PAS but terminates at Alu/unique boundary.
        # Requires: protein_coding intron + GN tag match + Alu boundary.
        if not pas:
            intron_hits_alu = parser.query_intron(chrom, ref_end, strand)
            if intron_hits_alu:
                exon_hits_alu = parser.query_exon(chrom, ref_end, strand)
                if not exon_hits_alu:
                    intron_gene_alu = intron_hits_alu[0][0]
                    if gn == intron_gene_alu:
                        is_alu, alu_label = check_alu_boundary(
                            fasta, chrom, ref_end, strand, alu_trees)
                        if is_alu:
                            site_id = registry.get_alu(
                                intron_gene_alu, chrom, ref_end, strand)
                            rows.append({
                                'cb':        row.cb,
                                'ub':        row.ub,
                                'gn':        gn,
                                'chrom':     chrom,
                                'ref_end':   ref_end,
                                'strand':    strand,
                                'site_id':   site_id,
                                'site_type': 'AluCE',
                                'pas_signal': alu_label,
                                'gene':      intron_gene_alu,
                            })
            continue  # no PAS and not Alu → skip

        # ── 4. Classify: CE or PA ─────────────────────────────────────────
        intron_hits = parser.query_intron(chrom, ref_end, strand)
        if intron_hits:
            exon_hits = parser.query_exon(chrom, ref_end, strand)
            if not exon_hits:
                intron_gene = intron_hits[0][0]
                if gn == intron_gene:
                    site_id   = registry.get_ce(intron_gene, chrom, ref_end, strand)
                    site_type = 'CE'
                    gene      = intron_gene
                else:
                    continue  # GN tag mismatch → skip
            else:
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
    parser.add_argument("--alu-bed",    default=None,
                        help="RepeatMasker Alu BED file for Alu CE detection "
                             "(chrom start end name score strand). "
                             "Generate with: gunzip -c rmsk.txt.gz | "
                             "awk '$13==\"Alu\"' | awk '{print $6\"\\t\"$7\"\\t\"$8\"\\t\"$11\"\\t\"$2\"\\t\"$10}' "
                             "| sort -k1,1 -k2,2n > alu_hg38.bed. "
                             "If not provided, falls back to hg38 soft-masking check.")
    parser.add_argument("--known-window", type=int, default=50,
                        help="Read scatter window for known 3' end matching (default: 50bp).")
    parser.add_argument("--batch-size", type=int, default=500_000)
    parser.add_argument("--include-no-ref", action="store_true",
                        help="Include NO_REF reads. Default: VALID_PAS only.")
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

    # Load Alu trees if provided
    alu_trees = None
    if args.alu_bed:
        alu_trees = load_alu_trees(args.alu_bed)
        print(f"[IsoCAPE] Alu CE detection: RepeatMasker BED ✅")
    else:
        print(f"[IsoCAPE] Alu CE detection: soft-masking fallback "
              f"(provide --alu-bed for higher precision)")

    valid_labels = {'VALID_PAS', 'NO_REF'} if args.include_no_ref else {'VALID_PAS'}

    total_reads = 0
    total_valid = 0
    total_sites = 0
    batch_num   = 0

    print(f"[IsoCAPE] Streaming (batch={args.batch_size:,}, known_window=±{args.known_window}bp)")
    print(f"[IsoCAPE] Output: known + CE + AluCE + PA")

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
            known_window=args.known_window,
            alu_trees=alu_trees
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
              f"CE={registry.n_ce:,} AluCE={registry.n_alu:,} PA={registry.n_pa:,} | "
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
    print(f"  Unique AluCE sites: {registry.n_alu:,}")
    print(f"  Unique PA sites:    {registry.n_pa:,}")
    print(f"  Total unique sites: {registry.n_unique:,}")
    print(f"  IsoCAPE rate:       {total_sites/total_valid:.3%} of VALID_PAS")
    print(f"  Output:             {args.out}")


if __name__ == "__main__":
    main()
