#!/usr/bin/env python3
"""
IsoCAPE: Build Count Matrix
-----------------------------
Takes annotated cryptic sites parquet from site_annotator.py,
performs positional clustering to merge nearby sites,
UMI deduplication per cell, and outputs a scanpy-ready AnnData.

Pipeline:
  1. Stream parquet → collect (gene, chrom, strand, ref_end) per read
  2. Positional clustering per gene (±tolerance bp)
  3. UMI deduplication per (cell, site)
  4. Build sparse cells × sites count matrix
  5. Add site metadata (site_type, pas_signal, chrom, rep_coord)
  6. Label lookup from cape_labels.csv
  7. Save AnnData (.h5ad)

Downstream:
    adata = sc.read_h5ad("IDC_cape.h5ad")
    sc.pp.filter_genes(adata, min_cells=3)
    ce = adata[:, adata.var['site_type'] == 'CE']

Usage:
python isocape/scripts/build_matrix.py \
    --parquet  results/IDC_cryptic.parquet \
    --out      results/IDC_cape.h5ad \
    [--labels  data/cape_labels.csv] \
    [--tolerance 10] \
    [--batch-size 500000]
"""

import argparse
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scipy.sparse as sp
import anndata as ad


# ---------------------------------------------------------------------------
# Positional clustering
# ---------------------------------------------------------------------------

def cluster_positions(positions_data, tolerance=10):
    """
    Cluster nearby ref_end positions within ±tolerance bp.
    Uses sliding window: compare each new position to the LAST
    position in current cluster (not the first) to avoid splitting
    dense clusters.

    Parameters
    ----------
    positions_data : list of (ref_end, site_id, pas_signal, ub, cb)
    tolerance      : int, bp window for merging

    Returns
    -------
    dict: rep_coord → list of (site_id, pas_signal, ub, cb)
    """
    if not positions_data:
        return {}

    sorted_data = sorted(positions_data, key=lambda x: x[0])
    clusters = {}
    current_items = [sorted_data[0]]

    for item in sorted_data[1:]:
        pos      = item[0]
        last_pos = current_items[-1][0]  # ← compare to LAST, not first

        if pos - last_pos <= tolerance:
            current_items.append(item)
        else:
            rep = int(np.median([x[0] for x in current_items]))
            clusters[rep] = current_items
            current_items = [item]

    rep = int(np.median([x[0] for x in current_items]))
    clusters[rep] = current_items
    return clusters


# ---------------------------------------------------------------------------
# PolyASite / PolyA_DB loader for PA site validation
# ---------------------------------------------------------------------------

def load_polyadb(bed_path, window=50):
    """
    Load PolyASite 2.0 / PolyA_DB bed file into an interval index.

    Chromosome normalization: adds 'chr' prefix if missing (1 → chr1, MT → chrM).

    Returns
    -------
    db : {(chrom, strand): sorted list of positions}
    """
    from collections import defaultdict

    db = defaultdict(list)
    n  = 0

    def norm_chrom(c):
        if not c.startswith('chr'):
            c = 'chrM' if c == 'MT' else f'chr{c}'
        return c

    with open(bed_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 6:
                continue
            chrom  = norm_chrom(fields[0])
            pos    = int(fields[2])  # use end as polyA site position
            strand = fields[5]
            db[(chrom, strand)].append(pos)
            n += 1

    for key in db:
        db[key].sort()

    print(f"[PolyADB] Loaded {n:,} polyA sites from {bed_path}")
    return db


def match_polyadb(chrom, pos, strand, db, window=50):
    """
    Check if a position is within ±window bp of any known polyA site.
    Uses binary search for O(log n).
    """
    import bisect
    sites = db.get((chrom, strand), [])
    if not sites:
        return False
    idx = bisect.bisect_left(sites, pos - window)
    while idx < len(sites) and sites[idx] <= pos + window:
        if abs(sites[idx] - pos) <= window:
            return True
        idx += 1
    return False




def collect_sites(parquet_path, batch_size=500_000):
    """
    Stream parquet and collect all reads grouped by
    (gene, chrom, strand, site_type).

    known/CE/PA are grouped separately to prevent cross-type
    clustering in Pass 2 (e.g. known reads merging with CE reads).

    Returns
    -------
    gene_data : {(gene, chrom, strand, site_type): [(ref_end, site_id, pas_signal, ub, cb), ...]}
    """
    pf        = pq.ParquetFile(parquet_path)
    gene_data = defaultdict(list)
    total     = 0

    print(f"[Pass 1] Collecting sites from {parquet_path} ...")

    for batch in pf.iter_batches(batch_size=batch_size):
        df = batch.to_pandas()
        total += len(df)

        for row in df.itertuples(index=False):
            # Include site_type in key → known/CE/PA never merge
            key = (row.gene, row.chrom, row.strand, row.site_type)
            gene_data[key].append((
                row.ref_end,
                row.site_id,
                row.pas_signal,
                row.ub,
                row.cb,
            ))

        if total % 5_000_000 == 0:
            print(f"  [Pass 1] {total:,} reads collected, "
                  f"{len(gene_data):,} gene-chrom-strand-type groups")

    print(f"[Pass 1] Done. {total:,} reads, {len(gene_data):,} groups")
    return gene_data


# ---------------------------------------------------------------------------
# Pass 2: cluster + UMI dedup → count matrix
# ---------------------------------------------------------------------------

def build_counts(gene_data, tolerance=10, min_reads=3, min_reads_alu=5):
    """
    For each (gene, chrom, strand, site_type) group:
    1. Cluster positions within ±tolerance bp
    2. Filter clusters with < min_reads total reads
    3. Assign coord-based site_id:
       CE     → GENE_CE_{rep_coord}
       AluCE  → GENE_AluCE_{rep_coord}
       PA     → GENE_PA_{rep_coord}
       known  → GENE_known_{rep_coord}
    4. UMI-deduplicate per (cb, site_id)
    5. Count UMIs per (cb, site_id)

    known/CE/PA are always clustered separately (enforced by key).
    """
    counts    = defaultdict(lambda: defaultdict(int))
    site_meta = {}

    print(f"[Pass 2] Clustering and counting "
          f"(tolerance={tolerance}bp, min_reads={min_reads}) ...")
    print(f"[Pass 2] known/CE/PA clustered separately (no cross-type merging)")

    n_groups      = len(gene_data)
    n_processed   = 0
    n_sites_total = 0
    n_filtered    = 0

    for (gene, chrom, strand, site_type), reads in gene_data.items():

        clusters = cluster_positions(reads, tolerance=tolerance)

        for rep_coord, items in clusters.items():

            threshold = min_reads_alu if site_type == 'AluCE' else min_reads
            if len(items) < threshold:
                n_filtered += 1
                continue

            pas_sigs   = [x[2] for x in items if x[2]]
            pas_signal = pas_sigs[0] if pas_sigs else None

            # Coord-based stable naming — site_type from key (never ambiguous)
            if site_type == 'CE':
                site_id = f"{gene}_CE_{rep_coord}"
            elif site_type == 'AluCE':
                site_id = f"{gene}_AluCE_{rep_coord}"
            elif site_type == 'PA':
                site_id = f"{gene}_PA_{rep_coord}"
            else:
                site_id = f"{gene}_known_{rep_coord}"

            if site_id not in site_meta:
                site_meta[site_id] = {
                    'gene':       gene,
                    'chrom':      chrom,
                    'strand':     strand,
                    'rep_coord':  rep_coord,
                    'site_type':  site_type,
                    'pas_signal': pas_signal,
                    'n_reads':    len(items),
                }
                n_sites_total += 1

            # UMI deduplication per cell
            seen_umis = defaultdict(set)
            for (ref_end, sid, pas, ub, cb) in items:
                if ub not in seen_umis[cb]:
                    seen_umis[cb].add(ub)
                    counts[cb][site_id] += 1

        n_processed += 1
        if n_processed % 10_000 == 0:
            print(f"  [Pass 2] {n_processed:,}/{n_groups:,} groups | "
                  f"{n_sites_total:,} sites kept | {n_filtered:,} filtered")

    print(f"[Pass 2] Done. {n_sites_total:,} sites kept, "
          f"{n_filtered:,} filtered (< {min_reads} reads)")
    return counts, site_meta


# ---------------------------------------------------------------------------
# Build AnnData
# ---------------------------------------------------------------------------

def build_anndata(counts, site_meta, labels_path=None, sample=None, polyadb=None):
    """
    Convert counts dict to sparse AnnData.

    obs  = cells (CB)
    var  = sites (site_id)
    X    = UMI count matrix (sparse)

    Parameters
    ----------
    sample  : str or None
        Sample name prefix for barcodes.
    polyadb : dict or None
        Output of load_polyadb(). If provided, PA sites are annotated as
        'PA_validated' (within ±50bp of known polyA site) or 'PA_novel'.
    """
    print(f"[Matrix] Building AnnData ...")

    # Normalize barcodes
    def normalize_bc(bc, sample):
        bc = bc.replace('-1', '').replace('-', '')
        if sample:
            return f"{sample}_{bc}"
        return bc

    all_cells_raw  = sorted(counts.keys())
    all_cells      = [normalize_bc(bc, sample) for bc in all_cells_raw]
    bc_map         = {raw: norm for raw, norm in zip(all_cells_raw, all_cells)}
    all_sites      = sorted(site_meta.keys())

    cell_idx = {bc: i for i, bc in enumerate(all_cells)}
    site_idx = {sid: j for j, sid in enumerate(all_sites)}

    rows, cols, vals = [], [], []
    for cb_raw, site_counts in counts.items():
        cb = bc_map[cb_raw]
        i  = cell_idx[cb]
        for sid, count in site_counts.items():
            if sid in site_idx:
                j = site_idx[sid]
                rows.append(i)
                cols.append(j)
                vals.append(count)

    X = sp.csr_matrix(
        (vals, (rows, cols)),
        shape=(len(all_cells), len(all_sites)),
        dtype=np.float32
    )

    # obs (cells)
    obs = pd.DataFrame(index=all_cells)
    obs.index.name = None

    # var (sites)
    var = pd.DataFrame(index=all_sites)
    var.index.name = None
    for col in ['gene', 'chrom', 'strand', 'rep_coord', 'site_type', 'pas_signal']:
        var[col] = [site_meta[sid].get(col) for sid in all_sites]
    var['feature_types'] = 'IsoCAPE'

    # PolyA_DB validation for PA sites
    if polyadb:
        print(f"[Matrix] Annotating PA sites against PolyA_DB ...")
        pa_mask = var['site_type'] == 'PA'
        n_validated = 0
        site_types_updated = var['site_type'].copy()
        for sid in var.index[pa_mask]:
            chrom = var.loc[sid, 'chrom']
            pos   = int(var.loc[sid, 'rep_coord'])
            strand = var.loc[sid, 'strand']
            if match_polyadb(chrom, pos, strand, polyadb, window=50):
                site_types_updated[sid] = 'PA_validated'
                n_validated += 1
            else:
                site_types_updated[sid] = 'PA_novel'
        var['site_type'] = site_types_updated
        n_novel = pa_mask.sum() - n_validated
        print(f"[Matrix] PA sites: {n_validated:,} validated, {n_novel:,} novel")

    # Label lookup
    if labels_path and os.path.exists(labels_path):
        labels = pd.read_csv(labels_path)
        label_map = labels.set_index('site_id')['label'].to_dict()
        clin_map  = labels.set_index('site_id')['clinical_relevance'].to_dict()
        var['label']               = var.index.map(label_map).fillna('unknown')
        var['clinical_relevance']  = var.index.map(clin_map).fillna('—')
        print(f"[Matrix] Applied labels to {var['label'].ne('unknown').sum()} sites")
    else:
        var['label']              = 'unknown'
        var['clinical_relevance'] = '—'

    adata = ad.AnnData(X=X, obs=obs, var=var)

    print(f"[Matrix] AnnData: {adata.n_obs:,} cells × {adata.n_vars:,} sites")
    print(f"  known sites:      {(var['site_type']=='known').sum():,}")
    print(f"  CE sites:         {(var['site_type']=='CE').sum():,}")
    print(f"  AluCE sites:      {(var['site_type']=='AluCE').sum():,}")
    if polyadb:
        print(f"  PA_validated:     {(var['site_type']=='PA_validated').sum():,}")
        print(f"  PA_novel:         {(var['site_type']=='PA_novel').sum():,}")
    else:
        print(f"  PA sites:         {(var['site_type']=='PA').sum():,}")
    print(f"  Total UMIs:       {int(X.sum()):,}")

    return adata


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="IsoCAPE: Build count matrix from cryptic sites parquet"
    )
    parser.add_argument("--parquet",    required=True,
                        help="Cryptic sites parquet from site_annotator.py")
    parser.add_argument("--out",        required=True,
                        help="Output AnnData h5ad path")
    parser.add_argument("--labels",     default=None,
                        help="cape_labels.csv for clinical annotation")
    parser.add_argument("--polyadb",    default=None,
                        help="PolyASite 2.0 / PolyA_DB bed file for PA site validation. "
                             "PA sites within ±50bp of a known polyA site → 'PA_validated'. "
                             "PA sites with no match → 'PA_novel'. "
                             "Recommended: ~/reference/hg38/polyasite2_hg38.bed")
    parser.add_argument("--sample",     default=None,
                        help="Sample name prefix for barcodes (e.g. 'IDC'). "
                             "Output barcodes: SAMPLE_BARCODE (strips -1 suffix). "
                             "Use same name as IsoDecipher for concat compatibility.")
    parser.add_argument("--tolerance",  type=int, default=10,
                        help="Positional clustering tolerance in bp (default: 10)")
    parser.add_argument("--min-reads",  type=int, default=3,
                        help="Minimum total reads per site cluster across all cells (default: 3). "
                             "Removes absolute singleton sites. "
                             "For cell-level filtering use sc.pp.filter_genes(adata, min_cells=3) downstream.")
    parser.add_argument("--min-reads-alu", type=int, default=5,
                        help="Min reads per AluCE cluster (default: 5, stricter than CE).")
    parser.add_argument("--batch-size", type=int, default=500_000,
                        help="Reads per batch for Pass 1 (default: 500,000)")
    return parser.parse_args()


def main():
    args = parse_args()

    # Pass 1: collect
    gene_data = collect_sites(args.parquet, batch_size=args.batch_size)

    # Pass 2: cluster + count
    counts, site_meta = build_counts(
        gene_data,
        tolerance=args.tolerance,
        min_reads=args.min_reads,
        min_reads_alu=args.min_reads_alu,
    )

    # Load PolyA_DB (optional)
    polyadb = None
    if args.polyadb:
        polyadb = load_polyadb(args.polyadb)

    # Build AnnData
    adata = build_anndata(
        counts, site_meta,
        labels_path=args.labels,
        sample=args.sample,
        polyadb=polyadb,
    )

    # Save
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    adata.write_h5ad(args.out)

    print(f"\n[SUCCESS] Saved → {args.out}")
    print(f"\nDownstream usage:")
    print(f"  import scanpy as sc")
    print(f"  adata = sc.read_h5ad('{args.out}')")
    print(f"  sc.pp.filter_genes(adata, min_cells=3)  # remove low-confidence sites")
    print(f"  ce = adata[:, adata.var['site_type'] == 'CE']  # cryptic exon sites only")


if __name__ == "__main__":
    main()
