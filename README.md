# IsoCAPE

**IsoCAPE** is a cryptic and premature polyadenylation scanner for RNA-seq data.

It detects polyadenylation sites that are **absent from GTF annotation** — intronic termination events, cryptic terminal exons, and novel 3' end extensions — from standard RNA-seq BAM files, with no specialized assay or library preparation required.

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue)]()
[![Compatible with Cell Ranger](https://img.shields.io/badge/input-any%20RNA--seq%20BAM-green)]()

---

## Concept

Standard RNA-seq pipelines assign reads to annotated gene boundaries. Polyadenylation events that occur inside introns — driven by cryptic terminal exons — are collapsed into total gene counts or discarded entirely. These events are not noise: they are the molecular mechanism behind some of the most clinically consequential isoforms in cancer.

IsoCAPE reinterprets the genomic termination positions of 3' RNA-seq reads to detect where transcripts end — including positions the annotation does not know about.

```
RNA-seq reads → Cell Ranger / STAR → gene counts (cryptic site lost)

Cell        AR
--------------------
cell_1      45
cell_2      62
cell_3      38
```

IsoCAPE resolves termination positions, compares them against the GTF, and identifies unannotated polyadenylation sites:

```
RNA-seq reads → Cell Ranger / STAR → gene counts (cryptic site lost)

Cell        AR
--------------------
cell_1      45
cell_2      62
cell_3      38
```

IsoCAPE detects reads terminating inside AR intron 3, verifies the upstream PAS signal, and reports the cryptic site:

```
Cell      AR_CE3 (AR-V7, cryptic intron 3 site)
------------------------------------------------
cell_1          1
cell_2         44
cell_3          1
```

cell_2 is expressing AR-V7 at high level. This is invisible to standard pipelines. Combine with IsoDecipher for the complete APA landscape:

```python
apa     = sc.read_h5ad("isodecipher_output.h5ad")  # AR_G0, AR_G1... (annotated)
cryptic = sc.read_h5ad("isocape_output.h5ad")       # AR_CE3 (cryptic)
combined = ad.concat([apa, cryptic], axis=1)
```

---

## The AR-V7 case

AR-V7 is the leading resistance mechanism to enzalutamide and abiraterone in castration-resistant prostate cancer (CRPC). Its molecular origin is premature polyadenylation driven by a cryptic terminal exon (CE3) in AR intron 3. CE3 carries its own polyadenylation signal (PAS), causing transcription to terminate mid-gene and produce a truncated, constitutively active AR variant that is insensitive to androgen receptor antagonists.

Current clinical detection of AR-V7 requires specialized liquid biopsy assays from CTCs, plasma, or urine — a separate test, a separate sample, a separate cost.

IsoCAPE detects the same signal from BAM files already sequenced. Every archived prostate cancer RNA-seq dataset is a retrospectively scannable source of AR-V7 status, at single-cell resolution.

```
AR intron 3
─────────────────────────────────────────────────
Exon 3          CE3 (cryptic)          Exon 4
────┤           ┌──────────────┐       ├────────
    │     ~~~~~~│ PAS: AATAAA  │~~~~~  │
    │           │ poly-A tail  │       │
    │           └──────────────┘       │
    │                ↑                 │
    │         3' reads terminate here  │
    │         → AR_CE3 pileup          │
    │         → AR-V7 isoform          │
─────────────────────────────────────────────────
```

CE3 is not in the standard GTF. Cell Ranger assigns these reads to AR total counts. IsoCAPE detects the intronic pileup, confirms the AATAAA signal upstream, and labels the site `AR_CE3`.

---

## Why IsoCAPE?

| Feature | scTail | CRYPTID-exon | Sierra | IsoCAPE |
|---------|--------|--------------|--------|---------|
| Works on standard reads 2 | ❌ | ✅ | ✅ | ✅ |
| Reads 1 required | ✅ (required) | ❌ | ❌ | ❌ |
| Single-cell resolution | ✅ | ❌ | ✅ | ✅ |
| Bulk RNA-seq | ❌ | ✅ | ✅ | ✅ |
| Internal priming exclusion | ✅ | N/A | ❌ | ✅ |
| Intronic / cryptic site detection | partial | ✅ | ❌ | ✅ |
| PAS signal verification (AATAAA + variants) | ❌ | ❌ | ❌ | ✅ |
| GTF-exclusion naming | ❌ | ❌ | ❌ | ✅ |
| Clinical label layer | ❌ | ❌ | ❌ | ✅ |
| Scanpy-ready output | ✅ | ❌ | ❌ | ✅ |
| Works on archived BAMs | ❌ | ✅ | ✅ | ✅ |

**Why reads 1 independence matters:** The majority of public 10x Genomics datasets have reads 1 trimmed or not sequenced. scTail, the only existing single-cell PAS detection tool, requires preserved reads 1 and cannot be applied to most archived datasets. IsoCAPE operates entirely on reads 2 — the same reads used by Cell Ranger — making it compatible with every standard 3' scRNA-seq dataset ever deposited.

**Why IsoCAPE's approach matters:** De novo tools like Sierra call any read pileup as a peak — producing coordinate-only names like `chr7:140833972` with no gene context and high noise. IsoCAPE uses the GTF as an exclusion reference to define what is already known (handled by IsoDecipher), then applies two-layer filtering — internal priming exclusion and PAS signal verification — to confirm only biologically meaningful unannotated sites. Every IsoCAPE site carries a gene name, a type (CE or PA), and a verified polyadenylation signal, making them immediately interpretable as clinical features or ML inputs without post-processing.

**Why the label layer matters:** Unknown sites are not useless. They are ML features with full rights, retrospectively labelable as biology catches up, and candidates for novel discovery. The label layer is a living CSV — `cape_labels.csv` — that maps site IDs to known clinical relevance. It ships with curated entries for AR-V7 and grows as the community contributes.

---

## How It Works

```
RNA-seq BAM (bulk / scRNA-seq / spatial)
              │
              ▼
       bam_to_parquet.py
       ├── CB + UB tag filter (single-cell)
       ├── Valid barcode filter (auto-detects TSV/CSV/gz format)
       ├── UMI deduplication
       ├── Phred QC on 3' window
       └── Internal priming exclusion
           (A-tract scan of downstream genomic sequence)
              │
              ▼  reads.parquet
              │
              ▼
       site_annotator.py
       ├── Skip: UNK_GENE (no gene context)
       ├── Skip: known 3' end ±10bp (IsoDecipher's territory)
       ├── Skip: no PAS signal (noise)
       ├── Intronic + PAS → CE site (GENE_CE1, GENE_CE2 ...)
       │   e.g. AR intron 3 → AR_CE3 (AR-V7)
       └── Genic + PAS    → PA site (GENE_PA1, GENE_PA2 ...)
               shorter UTR / longer UTR / alternative last exon
              │
              ▼  cryptic_sites.parquet
              │
              ▼
       build_matrix.py
       ├── cells × sites count matrix
       └── AnnData (.h5ad)
              │
              ▼
       Downstream (scanpy)
       ├── sc.pp.filter_genes(adata, min_cells=3)
       ├── Label lookup (cape_labels.csv)
       └── Integration with IsoDecipher output
```

**Why ±10bp window matches IsoDecipher:**
IsoCAPE uses a 10bp window to identify known 3' ends — identical to IsoDecipher's clustering tolerance. This ensures zero overlap: any site IsoDecipher would assign to a G-group is skipped by IsoCAPE, and vice versa.

**CE vs PA:**
- `CE` (Cryptic Exon): read terminates inside an annotated intron with PAS signal. The canonical case for premature polyadenylation — AR-V7's CE3 is detected here.
- `PA` (Alternative Polyadenylation): read terminates within a known gene but at an unannotated position — shorter or longer than the GTF-recorded 3' end, or in a non-last exon. All are genuine APA events invisible to IsoDecipher.

---

## Installation

```bash
mamba create -n isocape python=3.10 -y
mamba activate isocape
pip install -r requirements.txt
```

---

## Quick Start

### Step 1: Extract reads from BAM

```bash
python isocape/scripts/bam_to_parquet.py \
    --bam    data/patient_01.bam \
    --ref    data/hg38.fa \
    --out    results/patient_01_reads.parquet \
    --barcodes data/filtered_barcodes.csv \
    [--cores 4] \
    [--window 80] \
    [--min-phred 20] \
    [--priming-threshold 8]
```

Output columns: `cb`, `ub`, `gn`, `chrom`, `ref_end`, `strand`, `priming_label`

### Step 2: Annotate cryptic sites

```bash
python isocape/scripts/annotator/site_annotator.py \
    --parquet results/patient_01_reads.parquet \
    --gtf     data/Homo_sapiens.GRCh38.115.gtf \
    --db      data/Homo_sapiens.GRCh38.115.gtf.db \
    --ref     data/hg38.fa \
    --out     results/patient_01_cryptic.parquet \
    [--window 10]
```

Output columns: `cb`, `ub`, `gn`, `chrom`, `ref_end`, `strand`, `site_id`, `site_type`, `pas_signal`, `gene`

### Step 3: Build count matrix

```bash
python isocape/scripts/build_matrix.py \
    --sites  results/patient_01_cryptic.parquet \
    --labels isocape/cape_labels.csv \
    --out    results/patient_01_cape.h5ad
```

### Step 4: Downstream analysis (Python)

```python
import scanpy as sc
import anndata as ad

# IsoCAPE output
cape = sc.read_h5ad("results/patient_01_cape.h5ad")

# Filter low-confidence sites (same as scanpy standard)
sc.pp.filter_genes(cape, min_cells=3)

# CE sites only (intronic premature polyadenylation)
ce = cape[:, cape.var['site_type'] == 'CE']

# Check AR-V7 status per cell
ar_v7 = ce[:, ce.var_names.str.startswith('AR_CE')]

# Combine with IsoDecipher output
apa = sc.read_h5ad("isodecipher_output.h5ad")
combined = ad.concat([apa, cape], axis=1)
```

---

## Site naming convention

| Site type | Example ID | Meaning |
|-----------|-----------|---------|
| CE (Cryptic Exon) | `AR_CE3` | Intronic termination + PAS signal — premature polyadenylation |
| PA (Alternative PA) | `BRAF_PA1` | Genic termination at unannotated position + PAS signal — shorter UTR, longer UTR, or alternative last exon |

Both CE and PA sites are absent from GTF annotation. IsoDecipher handles annotated 3' ends (G-groups); IsoCAPE handles everything else within gene boundaries.

Coordinate information (`chrom`, `ref_end`, `strand`) is stored in `adata.var` for all sites, enabling IGV inspection and future annotation.

---

## Label layer

`cape_labels.csv` ships with curated entries for known clinically relevant cryptic sites:

```
site_id,    gene,  isoform,  label,                                    clinical_relevance
AR_CE3,     AR,    AR-V7,    Cryptic exon 3 premature polyadenylation,  Enzalutamide / abiraterone resistance in CRPC
AR_CE1,     AR,    AR-V1,    Cryptic exon 1,                            Androgen deprivation resistance; significance unclear
AR_CE2,     AR,    AR-V9,    Cryptic exon 2,                            Co-occurs with AR-V7; synergistic resistance
```

Unknown sites receive `label = unknown` and `clinical_relevance = —`. They are retained as first-class features in the count matrix.

The label file is a plain CSV — contributions welcome via pull request.

---

## Integration with IsoDecipher

IsoCAPE and IsoDecipher have a clean division of labour:

- **IsoDecipher** quantifies APA at annotated 3' ends — GTF-anchored, high-fidelity G-groups
- **IsoCAPE** detects APA at unannotated positions — CE and PA sites the GTF misses

The ±10bp window in IsoCAPE's site annotator matches IsoDecipher's clustering tolerance exactly, ensuring **zero overlap** between their outputs. Their AnnData objects share the same structure and concatenate directly:

```python
import anndata as ad
import scanpy as sc

apa  = sc.read_h5ad("isodecipher_output.h5ad")   # annotated APA (G-groups)
cape = sc.read_h5ad("isocape_output.h5ad")         # cryptic sites (CE + PA)

combined = ad.concat([apa, cape], axis=1)
# combined.var['feature_types'] distinguishes the two
```

IsoCAPE CE and PA sites are also candidate vocabulary entries for **IsoFormer** — the foundation model learns from both annotated and cryptic termination events across normal differentiation and malignancy.

---

## Repository Structure

```
IsoCAPE/
├── IsoCAPE/
│   └── scripts/
│       ├── bam_to_parquet.py           # Step 1: BAM → reads parquet (parallel streaming)
│       ├── bam_to_parquet_parallel.py  # Step 1: parallel version (multi-core)
│       └── annotator/
│           ├── gtf_parser.py           # GTF/DB index builder
│           └── site_annotator.py       # Step 2: cryptic site annotation (CE + PA)
├── notebooks/
│   └── 01_AR-V7_demo.ipynb            # Validation: AR-V7 detection in CRPC dataset
├── data/
│   └── cape_labels.csv                # Curated clinical label layer
├── results/
│   └── figures/
├── requirements.txt
└── README.md
```

---

## Future Directions

- **Mouse / multi-species support** — mm10 GTF + reference genome
- **Spatial transcriptomics** — Visium BAM compatibility (spot-level APA resolution across tissue sections)
- **IsoCAPE-DB** — a community database of validated cryptic polyadenylation sites with clinical annotations across cancer types
- **IsoFormer integration** — IsoCAPE sites form the novel-site vocabulary for BAM-level foundation model training

---

## License

MIT License  
© 2026 Rene Yu-Hong Cheng
