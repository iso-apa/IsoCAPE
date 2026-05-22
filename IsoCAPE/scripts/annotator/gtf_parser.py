#!/usr/bin/env python3
"""
IsoCAPE: GTF Parser
--------------------
Builds two lookup structures for site annotation:

1. known_3p_ends  — {chrom: [(pos, strand, gene_name, transcript_id), ...]}
2. intron_index   — {chrom: [(start, end, strand, gene_name, transcript_id), ...]}

Supports two input modes:
  - GTF file (.gtf or .gtf.gz)   → parsed directly
  - gffutils DB (.gtf.db)        → faster, recommended if available

Usage:
    from gtf_parser import GTFParser

    # With db (fast, recommended)
    parser = GTFParser(gtf_path="hg38.gtf", db_path="hg38.gtf.db")
    parser.build()

    # GTF only
    parser = GTFParser(gtf_path="hg38.gtf")
    parser.build()

    # Query
    matches = parser.query_3p_end("chr7", 140453136, "+", window=20)
    introns = parser.query_intron("chrX", 66765940, "-")
"""

import gzip
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class TranscriptModel:
    gene_name:     str
    transcript_id: str
    chrom:         str
    strand:        str
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


class GTFParser:
    """
    Parses GTF or gffutils DB and builds lookup structures for IsoCAPE.

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
        self.gene_types = gene_types or ["protein_coding"]

        self.known_3p_ends: Dict[str, List[Tuple[int, str, str, str]]] = defaultdict(list)
        self.intron_index:  Dict[str, List[Tuple[int, int, str, str, str]]] = defaultdict(list)
        self._transcripts:  Dict[str, TranscriptModel] = {}

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
        print(
            f"[GTFParser] Done. "
            f"{sum(len(v) for v in self.known_3p_ends.values()):,} known 3' ends | "
            f"{sum(len(v) for v in self.intron_index.values()):,} introns indexed"
        )

    def query_3p_end(self, chrom, pos, strand, window=20):
        matches = []
        for (end_pos, end_strand, gene_name, tx_id) in self.known_3p_ends.get(chrom, []):
            if end_strand == strand and abs(end_pos - pos) <= window:
                matches.append((gene_name, tx_id))
        return matches

    def query_intron(self, chrom, pos, strand):
        hits = []
        for (i_start, i_end, i_strand, gene_name, tx_id) in self.intron_index.get(chrom, []):
            if i_strand == strand and i_start <= pos < i_end:
                hits.append((gene_name, tx_id))
        return hits

    # ------------------------------------------------------------------
    # DB parsing (fast path)
    # ------------------------------------------------------------------

    def _parse_from_db(self) -> None:
        try:
            import gffutils
        except ImportError:
            print("[GTFParser] gffutils not installed — falling back to GTF parsing.")
            self._parse_gtf()
            return

        db = gffutils.FeatureDB(self.db_path)
        n_transcripts = 0

        for tx in db.features_of_type('transcript'):
            gene_type = (
                tx.attributes.get('gene_type',    [None])[0] or
                tx.attributes.get('gene_biotype', [None])[0] or ''
            )
            if self.gene_types and gene_type not in self.gene_types:
                continue

            gene_name     = (tx.attributes.get('gene_name', [None])[0] or
                             tx.attributes.get('gene_id',   ['UNK'])[0])
            transcript_id = tx.attributes.get('transcript_id', ['UNK'])[0]
            chrom         = tx.chrom
            strand        = tx.strand
            tx_key        = f"{chrom}::{transcript_id}"

            self._transcripts[tx_key] = TranscriptModel(
                gene_name=gene_name, transcript_id=transcript_id,
                chrom=chrom, strand=strand,
            )

            for exon in db.children(tx, featuretype='exon', order_by='start'):
                self._transcripts[tx_key].exons.append((exon.start - 1, exon.end))

            n_transcripts += 1
            if n_transcripts % 50_000 == 0:
                print(f"  [GTFParser] {n_transcripts:,} transcripts loaded ...")

        print(f"[GTFParser] Loaded {n_transcripts:,} transcripts from DB")

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

                attr_dict     = self._parse_attributes(attrs)
                gene_type     = attr_dict.get("gene_type") or attr_dict.get("gene_biotype", "")
                if self.gene_types and gene_type not in self.gene_types:
                    continue

                gene_name     = attr_dict.get("gene_name", attr_dict.get("gene_id", "UNK"))
                transcript_id = attr_dict.get("transcript_id", "UNK")
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
        for tx in self._transcripts.values():
            if not tx.exons:
                continue
            self.known_3p_ends[tx.chrom].append(
                (tx.three_prime_end(), tx.strand, tx.gene_name, tx.transcript_id)
            )
            for (i_start, i_end) in tx.introns():
                self.intron_index[tx.chrom].append(
                    (i_start, i_end, tx.strand, tx.gene_name, tx.transcript_id)
                )

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

    p = argparse.ArgumentParser(description="IsoCAPE GTF Parser — sanity check")
    p.add_argument("--gtf",    required=True)
    p.add_argument("--db",     default=None, help="gffutils DB path (optional, faster)")
    p.add_argument("--chrom",  default="chrX")
    p.add_argument("--pos",    type=int, default=66765940, help="AR locus default")
    p.add_argument("--strand", default="-")
    args = p.parse_args()

    gp = GTFParser(gtf_path=args.gtf, db_path=args.db)
    gp.build()

    print(f"\n--- Query: {args.chrom}:{args.pos} ({args.strand}) ---")
    matches = gp.query_3p_end(args.chrom, args.pos, args.strand)
    if matches:
        print(f"Known 3' end: {matches}")
    else:
        hits = gp.query_intron(args.chrom, args.pos, args.strand)
        if hits:
            print(f"Intronic (cryptic candidate): {hits}")
        else:
            print("Novel / intergenic")
