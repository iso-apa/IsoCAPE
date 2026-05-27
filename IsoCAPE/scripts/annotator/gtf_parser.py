#!/usr/bin/env python3
"""
IsoCAPE: GTF Parser (with interval tree for fast intron lookup)
---------------------------------------------------------------
Builds two lookup structures:

1. known_3p_ends  — {chrom: sorted list of (pos, strand, gene, tx_id)}
                    queried with binary search (bisect)

2. intron_trees   — {chrom: {strand: IntervalTree}}
                    queried with interval tree O(log n) instead of O(n)

Chromosome names are normalized to UCSC style (chr-prefix) so
Ensembl GTF ('1','X','MT') matches Cell Ranger BAM ('chr1','chrX','chrM').

Supports GTF file or gffutils DB input.

Usage:
    from gtf_parser import GTFParser

    parser = GTFParser(gtf_path="hg38.gtf", db_path="hg38.gtf.db")
    parser.build()

    matches = parser.query_3p_end("chr17", 41196312, "-", window=10)
    introns = parser.query_intron("chrX", 66765940, "-")
"""

import gzip
import bisect
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# TranscriptModel
# ---------------------------------------------------------------------------

@dataclass
class TranscriptModel:
    gene_name:     str
    transcript_id: str
    chrom:         str
    strand:        str
    tx_biotype:    str = 'protein_coding'  # transcript-level biotype
    exons:         List[Tuple[int, int]] = field(default_factory=list)

    def three_prime_end(self) -> int:
        if not self.exons:
            return -1
        if self.strand == '+':
            return max(e[1] for e in self.exons)
        else:
            return min(e[0] for e in self.exons)

    def introns(self) -> List[Tuple[int, int]]:
        if len(self.exons) < 2:
            return []
        sorted_exons = sorted(self.exons, key=lambda e: e[0])
        return [
            (sorted_exons[i][1], sorted_exons[i + 1][0])
            for i in range(len(sorted_exons) - 1)
        ]


# ---------------------------------------------------------------------------
# Simple interval tree (no external dependency)
# ---------------------------------------------------------------------------

class IntervalNode:
    """Node in an augmented interval tree (center-based)."""
    __slots__ = ['center', 'intervals', 'left', 'right']

    def __init__(self, center):
        self.center    = center
        self.intervals = []   # [(start, end, data), ...]
        self.left      = None
        self.right     = None


class IntervalTree:
    """
    Static interval tree for fast overlap queries.
    Build once, query many times.
    Each interval stores (start, end, data).
    query(pos) returns all intervals containing pos.
    """

    def __init__(self):
        self._intervals = []
        self._root      = None
        self._built     = False

    def add(self, start: int, end: int, data):
        self._intervals.append((start, end, data))

    def build(self):
        if self._intervals:
            self._root = self._build(self._intervals)
        self._built = True

    def _build(self, intervals):
        if not intervals:
            return None
        # Choose center as median of interval midpoints
        center = sorted((s + e) // 2 for s, e, _ in intervals)[len(intervals) // 2]
        node   = IntervalNode(center)
        left_ivs  = []
        right_ivs = []
        for iv in intervals:
            s, e, d = iv
            if e < center:
                left_ivs.append(iv)
            elif s > center:
                right_ivs.append(iv)
            else:
                node.intervals.append(iv)
        node.left  = self._build(left_ivs)
        node.right = self._build(right_ivs)
        return node

    def query(self, pos: int) -> list:
        results = []
        self._query(self._root, pos, results)
        return results

    def _query(self, node, pos, results):
        if node is None:
            return
        if pos < node.center:
            for s, e, d in node.intervals:
                if s <= pos:
                    results.append(d)
            self._query(node.left, pos, results)
        elif pos > node.center:
            for s, e, d in node.intervals:
                if e >= pos:  # ← include endpoint (closed interval)
                    results.append(d)
            self._query(node.right, pos, results)
        else:
            for s, e, d in node.intervals:
                results.append(d)

    def __len__(self):
        return len(self._intervals)


# ---------------------------------------------------------------------------
# GTFParser
# ---------------------------------------------------------------------------

class GTFParser:
    """
    Parses GTF or gffutils DB and builds fast lookup structures.

    Parameters
    ----------
    gtf_path   : str
    db_path    : str or None — gffutils .db file (faster)
    gene_types : list of str — biotypes to include. None = all.
    """

    def __init__(
        self,
        gtf_path:   str,
        db_path:    Optional[str] = None,
        gene_types: Optional[List[str]] = None,
    ):
        self.gtf_path   = gtf_path
        self.db_path    = db_path

        # gene_types for known_3p_ends and intron_trees
        # Strict: only protein_coding transcripts define known APA sites and CE candidates
        self.gene_types = gene_types or ["protein_coding"]

        # exon_types for exon_trees — inclusive but conservative
        # Purpose: prevent false CE calls
        # Include: protein_coding related biotypes only
        # Exclude: lncRNA (too aggressive, overlaps protein_coding introns)
        #          retained_intron (IS the intron — would kill real CE calls)
        #          TEC (unconfirmed)
        self.exon_types = {
            "protein_coding",
            "protein_coding_CDS_not_defined",
            "protein_coding_LoF",
            "nonsense_mediated_decay",  # NMD isoforms of protein_coding genes
            "non_stop_decay",           # NSD isoforms of protein_coding genes
        }

        # known_3p_ends: {chrom: sorted list of (pos, strand, gene, tx_id)}
        # sorted by pos for binary search
        self.known_3p_ends: Dict[str, List] = defaultdict(list)

        # intron_trees: {chrom: {strand: IntervalTree}}
        self.intron_trees: Dict[str, Dict[str, IntervalTree]] = defaultdict(
            lambda: {'+': IntervalTree(), '-': IntervalTree()}
        )

        # exon_trees: {chrom: {strand: IntervalTree}}
        # Used to check if a position is in ANY exon across all transcripts
        # CE requires: intronic + NOT in any exon
        self.exon_trees: Dict[str, Dict[str, IntervalTree]] = defaultdict(
            lambda: {'+': IntervalTree(), '-': IntervalTree()}
        )

        self._transcripts: Dict[str, TranscriptModel] = {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def build(self) -> None:
        if self.db_path:
            print(f"[GTFParser] Loading from DB: {self.db_path}")
            self._parse_from_db()
        else:
            print(f"[GTFParser] Parsing GTF: {self.gtf_path}")
            self._parse_gtf()
        self._build_indexes()
        n_ends    = sum(len(v) for v in self.known_3p_ends.values())
        n_introns = sum(
            len(tree)
            for strands in self.intron_trees.values()
            for tree in strands.values()
        )
        print(f"[GTFParser] Done. {n_ends:,} known 3' ends | {n_introns:,} introns indexed")

    def query_3p_end(self, chrom, pos, strand, window=10):
        """
        Binary search for known 3' ends within ±window bp.
        Returns list of (gene_name, transcript_id).
        """
        ends = self.known_3p_ends.get(chrom, [])
        if not ends:
            return []

        # Binary search for window
        lo = bisect.bisect_left(ends,  (pos - window,))
        hi = bisect.bisect_right(ends, (pos + window,))

        matches = []
        for i in range(lo, hi):
            end_pos, end_strand, gene, tx_id = ends[i]
            if end_strand == strand and abs(end_pos - pos) <= window:
                matches.append((gene, tx_id))
        return matches

    def query_intron(self, chrom, pos, strand):
        """
        Interval tree lookup for introns containing pos.
        Returns list of (gene_name, transcript_id).
        O(log n + k) where k = number of results.
        """
        trees = self.intron_trees.get(chrom)
        if not trees:
            return []
        tree = trees.get(strand)
        if not tree or not tree._built:
            return []
        return tree.query(pos)

    def query_exon(self, chrom, pos, strand):
        """
        Interval tree lookup for exons containing pos.
        Returns list of (gene_name, transcript_id).
        Used to verify CE sites: a position must be intronic AND
        NOT in any exon across all transcripts.
        """
        trees = self.exon_trees.get(chrom)
        if not trees:
            return []
        tree = trees.get(strand)
        if not tree or not tree._built:
            return []
        return tree.query(pos)

    # ------------------------------------------------------------------
    # Chrom normalization
    # ------------------------------------------------------------------

    @staticmethod
    def _normalize_chrom(chrom: str) -> str:
        """Ensembl → UCSC: '1' → 'chr1', 'X' → 'chrX', 'MT' → 'chrM'."""
        if chrom.startswith('chr'):
            return chrom
        if chrom == 'MT':
            return 'chrM'
        return f"chr{chrom}"

    # ------------------------------------------------------------------
    # DB parsing (fast)
    # ------------------------------------------------------------------

    def _parse_from_db(self) -> None:
        try:
            import gffutils
        except ImportError:
            print("[GTFParser] gffutils not installed — falling back to GTF.")
            self._parse_gtf()
            return

        db = gffutils.FeatureDB(self.db_path)
        n  = 0
        n_exon_only = 0

        for tx in db.features_of_type('transcript'):
            # transcript_biotype for exon_types check (transcript-level)
            tx_biotype = (
                tx.attributes.get('transcript_biotype', [None])[0] or
                tx.attributes.get('gene_type',          [None])[0] or
                tx.attributes.get('gene_biotype',       [None])[0] or ''
            )
            # gene_biotype for strict filter (gene-level protein_coding)
            gene_biotype = (
                tx.attributes.get('gene_biotype', [None])[0] or
                tx.attributes.get('gene_type',    [None])[0] or ''
            )

            gene_name     = (tx.attributes.get('gene_name', [None])[0] or
                             tx.attributes.get('gene_id',   ['UNK'])[0])
            transcript_id = tx.attributes.get('transcript_id', ['UNK'])[0]
            chrom         = self._normalize_chrom(tx.chrom)
            strand        = tx.strand
            tx_key        = f"{chrom}::{transcript_id}"

            # exon_trees: use transcript_biotype (inclusive)
            if tx_biotype in self.exon_types:
                for exon in db.children(tx, featuretype='exon', order_by='start'):
                    self.exon_trees[chrom][strand].add(
                        exon.start - 1, exon.end,
                        (gene_name, transcript_id)
                    )
                n_exon_only += 1

            # known_3p_ends + intron_trees: use gene_biotype (strict)
            if self.gene_types and gene_biotype not in self.gene_types:
                continue

            self._transcripts[tx_key] = TranscriptModel(
                gene_name=gene_name, transcript_id=transcript_id,
                chrom=chrom, strand=strand,
                tx_biotype=tx_biotype,
            )
            for exon in db.children(tx, featuretype='exon', order_by='start'):
                self._transcripts[tx_key].exons.append((exon.start - 1, exon.end))

            n += 1
            if n % 50_000 == 0:
                print(f"  [GTFParser] {n:,} transcripts loaded ...")

        print(f"[GTFParser] Loaded {n:,} protein_coding transcripts from DB")
        print(f"[GTFParser] Exon index covers {n_exon_only:,} transcripts "
              f"({', '.join(sorted(self.exon_types))})")

    # ------------------------------------------------------------------
    # GTF parsing (fallback)
    # ------------------------------------------------------------------

    def _parse_gtf(self) -> None:
        opener = gzip.open if self.gtf_path.endswith(".gz") else open
        n_exons = 0

        with opener(self.gtf_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                chrom, _, feature, start, end, _, strand, _, attrs = fields
                if feature != "exon":
                    continue

                attr_dict = self._parse_attributes(attrs)
                gene_type = attr_dict.get("gene_type") or attr_dict.get("gene_biotype", "")
                if self.gene_types and gene_type not in self.gene_types:
                    continue

                gene_name     = attr_dict.get("gene_name", attr_dict.get("gene_id", "UNK"))
                transcript_id = attr_dict.get("transcript_id", "UNK")
                chrom         = self._normalize_chrom(chrom)
                tx_key        = f"{chrom}::{transcript_id}"

                if tx_key not in self._transcripts:
                    self._transcripts[tx_key] = TranscriptModel(
                        gene_name=gene_name, transcript_id=transcript_id,
                        chrom=chrom, strand=strand,
                    )
                self._transcripts[tx_key].exons.append((int(start) - 1, int(end)))
                n_exons += 1

        print(f"[GTFParser] Loaded {len(self._transcripts):,} transcripts, {n_exons:,} exons")

    # ------------------------------------------------------------------
    # Index builder
    # ------------------------------------------------------------------

    def _build_indexes(self) -> None:
        print(f"[GTFParser] Building indexes ...")
        for tx in self._transcripts.values():
            if not tx.exons:
                continue
            chrom = tx.chrom

            # known 3' ends — only protein_coding transcripts
            # Exclude retained_intron, NMD etc. from known 3'end index
            # to prevent intron reads being mis-labeled as 'known'
            if tx.tx_biotype == 'protein_coding':
                self.known_3p_ends[chrom].append(
                    (tx.three_prime_end(), tx.strand, tx.gene_name, tx.transcript_id)
                )

            # intron interval trees (protein_coding only)
            for (i_start, i_end) in tx.introns():
                self.intron_trees[chrom][tx.strand].add(
                    i_start, i_end, (tx.gene_name, tx.transcript_id)
                )

        # Sort known_3p_ends by position for binary search
        for chrom in self.known_3p_ends:
            self.known_3p_ends[chrom].sort(key=lambda x: x[0])

        # Build all interval trees
        print(f"[GTFParser] Building interval trees ...")
        for chrom, strands in self.intron_trees.items():
            for strand, tree in strands.items():
                tree.build()

        for chrom, strands in self.exon_trees.items():
            for strand, tree in strands.items():
                tree.build()

        n_introns = sum(
            len(t._intervals)
            for strands in self.intron_trees.values()
            for t in strands.values()
        )
        n_exons = sum(
            len(t._intervals)
            for strands in self.exon_trees.values()
            for t in strands.values()
        )
        print(f"[GTFParser] Done. {len(self.known_3p_ends):,} chroms | "
              f"{n_introns:,} introns indexed | {n_exons:,} exons indexed")

    @staticmethod
    def _parse_attributes(attr_string: str) -> dict:
        attrs = {}
        for item in attr_string.strip().split(";"):
            item = item.strip()
            if not item:
                continue
            if " " in item:
                key, _, value = item.partition(" ")
                attrs[key] = value.strip('"')
            elif "=" in item:
                key, _, value = item.partition("=")
                attrs[key] = value.strip('"')
        return attrs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    import time

    p = argparse.ArgumentParser(description="IsoCAPE GTF Parser — sanity check")
    p.add_argument("--gtf",    required=True)
    p.add_argument("--db",     default=None)
    p.add_argument("--chrom",  default="chrX")
    p.add_argument("--pos",    type=int, default=66765940)
    p.add_argument("--strand", default="-")
    args = p.parse_args()

    t0 = time.time()
    gp = GTFParser(gtf_path=args.gtf, db_path=args.db)
    gp.build()
    print(f"Build time: {time.time()-t0:.1f}s")

    t1 = time.time()
    for _ in range(10000):
        gp.query_intron(args.chrom, args.pos, args.strand)
    print(f"10,000 intron queries: {time.time()-t1:.3f}s")

    print(f"\n--- Query: {args.chrom}:{args.pos} ({args.strand}) ---")
    matches = gp.query_3p_end(args.chrom, args.pos, args.strand)
    if matches:
        print(f"Known 3' end: {matches}")
    else:
        hits = gp.query_intron(args.chrom, args.pos, args.strand)
        if hits:
            print(f"Intronic: {hits[:3]}")
        else:
            print("Novel / intergenic")
