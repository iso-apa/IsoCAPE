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

IsoCAPE resolves termination positions, compares them against the GTF, and labels intronic pileups:

```
Cell      AR_FL (annotated)   AR_CE3 (AR-V7)
----------------------------------------------
cell_1          44                  1
cell_2          18                 44
cell_3          37                  1
```

cell_2 is expressing AR-V7 at high level. This is invisible to standard pipelines. It is directly detectable from the existing BAM file.

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
| PAS signal verification (AATAAA) | ❌ | ❌ | ❌ | ✅ |
| GTF-hybrid naming | ❌ | ❌ | ❌ | ✅ |
| Clinical label layer | ❌ | ❌ | ❌ | ✅ |
| Scanpy-ready output | ✅ | ❌ | ❌ | ✅ |
| Works on archived BAMs | ❌ | ✅ | ✅ | ✅ |

**Why reads 1 independence matters:** The majority of public 10x Genomics datasets have reads 1 trimmed or not sequenced. scTail, the only existing single-cell PAS detection tool, requires preserved reads 1 and cannot be applied to most archived datasets. IsoCAPE operates entirely on reads 2 — the same reads used by Cell Ranger — making it compatible with every standard 3' scRNA-seq dataset ever deposited.

**Why the GTF-hybrid approach matters:** Fully de novo tools fragment reads across low-confidence peaks and produce coordinate-only feature names with no biological interpretability. Pure GTF tools miss everything outside annotation. IsoCAPE anchors known sites to GTF structure (preserving IsoDecipher compatibility) and assigns systematic, human-readable names to novel sites — `AR_CE3`, `GENE_novel_1` — that are immediately usable as ML features or clinical labels.

**Why the label layer matters:** Unknown sites are not useless. They are ML features with full rights, retrospectively labelable as biology catches up, and candidates for novel discovery. The label layer is a living CSV — `cape_labels.csv` — that maps site IDs to known clinical relevance. It ships with curated entries for AR-V7 and grows as the community contributes.

---

## How It Works

```
RNA-seq BAM (bulk / scRNA-seq / spatial)
              │
              ▼
       IsoCAPE Scanner
       ├── CB + UB tag filter (single-cell) or read-level (bulk)
       ├── UMI deduplication
       ├── Phred QC on 3' window
       └── Internal priming exclusion
           (A-tract scan of downstream genomic sequence via hg38)
              │
              ▼
       Site Annotation
       ├── GTF comparison (±20bp)
       │   ├── Matches annotated 3' end → IsoDecipher G-group (passthrough)
       │   └── Intronic / unannotated  → positional clustering
       │                                  ├── PAS signal check (AATAAA + 12 variants)
       │                                  ├── Minimum cell / UMI confidence filter
       │                                  └── Named: GENE_CE1, GENE_novel_1, ...
              │
              ▼
       Label Lookup (cape_labels.csv)
       ├── AR_CE3      → "AR-V7; enzalutamide/abiraterone resistance (CRPC)"
       ├── AR_CE1      → "AR-V1; unknown clinical significance"
       └── GENE_novel  → "unknown"
              │
              ▼
       Count Matrix (cells × sites)  +  AnnData (.h5ad)
```

---

## Installation

```bash
mamba create -n isocape python=3.10 -y
mamba activate isocape
pip install -r requirements.txt
```

---

## Quick Start

### Step 1: Scan BAM for cryptic termination sites

```bash
python isocape/scripts/scan_bam.py \
    --bam data/patient_01.bam \
    --ref data/hg38.fasta \
    --gtf data/Homo_sapiens.GRCh38.115.gtf \
    --out results/patient_01_sites.parquet \
    [--min-phred 20] \
    [--min-cells 10] \
    [--priming-threshold 8]
```

Output columns: `cb`, `ub`, `chrom`, `ref_end`, `strand`, `site_id`, `site_type`, `pas_signal`, `priming_label`

### Step 2: Build count matrix

```bash
python isocape/scripts/build_matrix.py \
    --sites results/patient_01_sites.parquet \
    --labels isocape/cape_labels.csv \
    --out results/patient_01_cape.h5ad
```

### Step 3: Downstream analysis (Python)

```python
import scanpy as sc

adata = sc.read_h5ad("results/patient_01_cape.h5ad")

# Filter to high-confidence novel sites
adata = adata[:, adata.var['site_type'] == 'novel']

# Check AR-V7 status per cell
ar_v7 = adata[:, adata.var['label'].str.contains('AR-V7', na=False)]
```

---

## Site naming convention

| Site type | Example ID | Meaning |
|-----------|-----------|---------|
| Annotated (GTF match) | `AR_G0` | Known 3'UTR end, IsoDecipher-compatible |
| Cryptic exon (labeled) | `AR_CE3` | Known cryptic exon with clinical label |
| Cryptic exon (novel) | `KRAS_CE1` | Intronic pileup + PAS signal, not in GTF |
| Novel 3'UTR extension | `CD59_novel_1` | Beyond annotated 3' end |

Coordinate information is stored in `adata.var` for all novel sites, enabling IGV inspection and future annotation.

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

IsoCAPE and IsoDecipher are designed to complement each other:

- **IsoDecipher** quantifies APA within annotated 3'UTR space (the known world)
- **IsoCAPE** detects polyadenylation outside annotation (what the GTF misses)

Their outputs share the same AnnData structure and can be concatenated along the `var` axis for a complete isoform landscape:

```python
import anndata as ad

gex    = sc.read_h5ad("isodecipher_output.h5ad")   # annotated APA
cryptic = sc.read_h5ad("isocape_output.h5ad")       # cryptic sites

combined = ad.concat([gex, cryptic], axis=1)
```

IsoCAPE novel sites are also the source vocabulary for **IsoDecipher-GPT** — the foundation model learns from both annotated and cryptic termination events.

---

## Repository Structure

```
IsoCAPE/
├── isocape/
│   ├── scripts/
│   │   ├── scan_bam.py          # Step 1: BAM → cryptic site detection
│   │   ├── annotate_sites.py    # Step 2: GTF comparison + site naming
│   │   ├── build_matrix.py      # Step 3: cells × sites count matrix
│   │   └── utils.py             # PAS signal detection, A-tract filter
│   └── __init__.py
├── notebooks/
│   └── 01_AR-V7_demo.ipynb      # Validation: AR-V7 detection in CRPC dataset
├── data/
│   └── cape_labels.csv          # Curated clinical label layer
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
