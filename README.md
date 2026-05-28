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

IsoCAPE detects reads terminating inside AR intron 3, verifies the upstream PAS signal, and reports the cryptic site:

```
Cell      AR_CE3 (AR-V7, cryptic intron 3 site)
------------------------------------------------
cell_1          1
cell_2         44
cell_3          1
```


cell_2 is expressing AR-V7 at high level. This is invisible to standard pipelines. IsoCAPE and IsoDecipher outputs are stored as separate `obsm` layers in the GEX AnnData object — see [Integration with IsoDecipher](#integration-with-isodecipher) for the complete integration pattern.


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
           (A-tract scan ≥8 consecutive A's in downstream 20bp;
            reads failing this check are labelled INTERNAL_PRIME
            and excluded from downstream steps)
              │
              ▼  reads.parquet  [priming_label = VALID_PAS | INTERNAL_PRIME | NO_REF]
              │
              ▼
       site_annotator.py
       ├── Skip: UNK_GENE (no gene context from GN tag)
       ├── Known 3' end check: ±50bp window against protein_coding
       │   transcript 3' ends (tx_biotype = protein_coding only;
       │   retained_intron / NMD transcript ends excluded)
       │   → match: labelled `known`; always output as reference signal
       ├── PAS signal check: scan 60bp upstream of ref_end
       │   PAS must be ≥10bp from ref_end (PAS_MIN_DIST)
       │   → no PAS: skip (noise)
       ├── CE classification (three independent filters required):
       │   1. query_intron(): read falls in protein_coding intron
       │   2. query_exon():   NOT in any exon of protein_coding /
       │                      NMD / non_stop_decay transcript
       │   3. GN tag == intron gene (neighboring-gene noise filter)
       │   → all three pass: labelled `CE`  (GENE_CE_{coord})
       └── PA classification (remaining genic + PAS reads):
           → labelled `PA`  (GENE_PA_{coord})
              │
              ▼  cryptic.parquet
              │
              ▼
       build_matrix.py
       ├── Pass 1: cluster reads per (gene, chrom, strand, site_type)
       │   within ±10bp tolerance; CE / PA / known never merge
       ├── Pass 2: filter clusters with < min_reads (default 3)
       ├── Coord-based stable naming: GENE_CE_{rep_coord}
       ├── PolyASite 2.0 validation (--polyadb):
       │   PA within ±50bp of known polyA site → PA_validated
       │   PA with no database match           → PA_novel
       └── AnnData (.h5ad): cells × sites
              │
              ▼
       Downstream (scanpy)
       ├── sc.pp.filter_genes(adata, min_cells=3)
       ├── CE fraction: CE / (CE + PA + known)  per gene
       ├── Label lookup (cape_labels.csv)
       └── Integration with IsoDecipher output (obsm layers)
```

**GTF biotype handling:**

| Index | Biotypes included | Purpose |
|-------|------------------|---------|
| `intron_trees` | `protein_coding` only | CE positive signal — intron must belong to protein-coding transcript |
| `exon_trees` | `protein_coding`, `protein_coding_CDS_not_defined`, `protein_coding_LoF`, `nonsense_mediated_decay`, `non_stop_decay` | CE negative filter — excludes positions exonic in any of these transcripts |
| `known_3p_ends` | `protein_coding` **transcript** (`tx_biotype`) | Known 3' end reference — retained_intron / NMD transcript ends excluded to prevent false known calls in intron regions |

**Known sites:**
`known` sites are always output alongside CE and PA. They represent reads terminating within ±50bp of a GTF-annotated protein-coding transcript 3' end. Use them as the denominator for CE fraction:

```
CE fraction = CE / (CE + PA + known)
```

**CE vs PA reliability:**

CE sites pass three independent filters (intron membership, exon exclusion, GN tag consistency) before PAS verification. PA sites pass only genic membership and PAS verification. **CE sites are substantially more reliable than PA sites.** See [Caveats](#caveats) for details.

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
    [--known-window 50]    # bp window for known 3' end matching (default: 50)
```

Output columns: `cb`, `ub`, `gn`, `chrom`, `ref_end`, `strand`, `site_id`, `site_type`, `pas_signal`, `gene`

`site_type` values: `CE` | `PA` | `known`

### Step 3: Build count matrix

```bash
python isocape/scripts/build_matrix.py \
    --parquet  results/patient_01_cryptic.parquet \
    --out      results/patient_01_cape.h5ad \
    --sample   patient_01 \               # sample prefix for barcodes (e.g. IDC)
    [--labels  isocape/cape_labels.csv] \
    [--polyadb data/polyasite2_hg38.bed] \ # PolyASite 2.0 for PA validation
    [--min-reads 3] \                      # minimum reads per site cluster
    [--tolerance 10]                       # clustering window in bp
```

`--sample`: formats barcodes as `SAMPLE_BARCODE` — use the same name as IsoDecipher for direct integration. Omit for standalone analysis.

`--polyadb`: when provided, PA sites are split into `PA_validated` (within ±50bp of a PolyASite 2.0 entry) and `PA_novel`. Download PolyASite 2.0 (hg38):

```bash
curl -L https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz \
     -o data/polyasite2_hg38.bed.gz && gunzip data/polyasite2_hg38.bed.gz
```

### Step 4: Downstream analysis (Python)

```python
import scanpy as sc
import numpy as np

# IsoCAPE output
cape = sc.read_h5ad("results/patient_01_cape.h5ad")

# Filter low-confidence sites
sc.pp.filter_genes(cape, min_cells=3)

# Site type breakdown
print(cape.var['site_type'].value_counts())
# CE              7,823
# PA_validated   31,241
# PA_novel       76,580
# known          37,512

# CE sites only (most reliable)
ce = cape[:, cape.var['site_type'] == 'CE']

# CE fraction per gene: CE / (CE + PA + known)
def ce_fraction(gene, cape):
    ce  = cape[:, (cape.var['gene']==gene) & (cape.var['site_type']=='CE')].X.sum()
    pa  = cape[:, (cape.var['gene']==gene) & cape.var['site_type'].isin(['PA','PA_validated','PA_novel'])].X.sum()
    kn  = cape[:, (cape.var['gene']==gene) & (cape.var['site_type']=='known')].X.sum()
    return float(ce / (ce + pa + kn)) if (ce + pa + kn) > 0 else 0

# Check AR-V7 status per cell
ar_v7 = ce[:, ce.var_names.str.startswith('AR_CE')]

# Integration with IsoDecipher (via obsm layers — no concat needed)
adata_gex = sc.read_h5ad("gex_annotated.h5ad")
# adata_gex.obsm['isocape']  → IsoCAPE matrix (cells × sites)
# adata_gex.obsm['isoform']  → IsoDecipher matrix (cells × G-groups)
# adata_gex.uns['isocape_features'] → site names
# adata_gex.uns['isocape_var']      → site metadata
```

---

## Site naming convention

| Site type | Example ID | Meaning |
|-----------|-----------|---------|
| `CE` | `MXD4_CE_2254615` | Intronic termination + PAS signal — premature polyadenylation. Three independent filters: intron membership, exon exclusion, GN tag consistency. Most reliable site type. |
| `PA_validated` | `ESR1_PA_152103175` | Genic + PAS, within ±50bp of a PolyASite 2.0 entry. Known alternative polyA site. |
| `PA_novel` | `TP53_PA_7676520` | Genic + PAS, no database match. Candidate novel APA site; treat with caution (see Caveats). |
| `known` | `FUS_known_31191575` | Read terminates within ±50bp of a GTF protein-coding transcript 3' end. Used as reference signal and CE fraction denominator. |

Site names embed the representative cleavage coordinate for cross-run stability. Coordinate-based naming ensures the same site receives the same ID across samples and pipeline versions.

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
- **IsoCAPE** detects APA at unannotated positions — CE and PA sites the GTF misses, plus `known` sites as reference signal

The recommended integration stores both modalities as `obsm` layers in the GEX AnnData object, avoiding axis-1 concatenation:

```python
import scanpy as sc
import scipy.sparse as sp
import numpy as np

# Load
adata_gex = sc.read_h5ad("gex_annotated.h5ad")   # GEX (preprocessed)
cape       = sc.read_h5ad("isocape_output.h5ad")   # IsoCAPE
apa        = sc.read_h5ad("isodecipher_output.h5ad") # IsoDecipher

# Align barcodes (IsoCAPE/IsoDecipher use SAMPLE_BARCODE format)
adata_gex.obs_names = ['SAMPLE_' + bc.replace('-1','') for bc in adata_gex.obs_names]

common_cape = adata_gex.obs_names.intersection(cape.obs_names)
common_apa  = adata_gex.obs_names.intersection(apa.obs_names)

# IsoCAPE layer
cell_to_idx = {c: i for i, c in enumerate(adata_gex.obs_names)}
cape_mat = sp.lil_matrix((adata_gex.n_obs, cape.n_vars), dtype='float32')
for cell in common_cape:
    cape_mat[cell_to_idx[cell]] = cape[cell].X
adata_gex.obsm['isocape'] = sp.csr_matrix(cape_mat)
adata_gex.uns['isocape_features'] = cape.var_names.tolist()
adata_gex.uns['isocape_var']      = cape.var.to_dict()

# IsoDecipher layer
iso_feats = apa.var_names[apa.var['feature_types'] == 'Isoform']
iso_mat = sp.lil_matrix((adata_gex.n_obs, len(iso_feats)), dtype='float32')
for cell in common_apa:
    iso_mat[cell_to_idx[cell]] = apa[cell, iso_feats].X
adata_gex.obsm['isoform'] = sp.csr_matrix(iso_mat)
adata_gex.uns['isoform_features'] = iso_feats.tolist()
```

IsoCAPE CE and PA sites are also candidate vocabulary entries for **IsoFormer** — the foundation model learns from both annotated and cryptic termination events across normal differentiation and malignancy.

---

## Caveats

### CE sites are more reliable than PA sites

CE sites pass **three independent filters** before PAS verification:

1. Must fall inside a `protein_coding` intron (`intron_trees`)
2. Must NOT overlap any exon of protein-coding, NMD, or non_stop_decay transcripts (`exon_trees`)
3. GN tag must match the intron's gene — neighboring-gene reads are excluded

PA sites require only: genic position + PAS signal + distance > 50bp from any known 3' end.

**Recommendation:** Use CE sites as the primary feature of interest. Apply stricter filtering for PA sites:

```python
# PA analysis: use only validated sites
pa_val = cape[:, cape.var['site_type'] == 'PA_validated']
sc.pp.filter_genes(pa_val, min_cells=10)  # stricter threshold
```

### Internal priming detection is incomplete

`bam_to_parquet` flags reads with ≥8 consecutive A's in the 20bp downstream genomic sequence as `INTERNAL_PRIME`. This catches obvious genomic polyA tracts but **does not capture**:

- Non-consecutive A-rich regions (e.g. `AATAAA…AAAAA` patterns)
- PAS signals that coincide with moderate A-rich context

Consequence: some CE and PA sites may reflect oligo-dT mispriming rather than genuine polyadenylation, particularly in A-rich intronic regions. The correlation between CE reads and GEX expression is a useful but imperfect diagnostic — genuine internal priming typically shows high CE/GEX correlation, but expression-linked CE events also exist.

**Future improvement:** sequence-context priming probability model (planned for IsoCAPE v2).

### Non-tumor reference is required for cancer-specificity assessment

IsoCAPE detects **CE usage**, not cancer-specific CE upregulation. Many CE sites are constitutive (present in all cell types) or reflect cell-type-specific APA rather than cancer biology. To identify cancer-specific events:

1. Compare CE fraction (CE / CE+PA+known) between tumor and non-tumor cells in the same dataset
2. Validate with IGV: cancer-specific sites should show a pileup in tumor BAM but not in matched normal or PBMC BAM
3. Use non-tumor cells within the dataset (macrophages, stromal cells) as an internal reference

Sites confirmed by both approaches (statistical enrichment + IGV) are high-confidence cancer-specific CE events.

### 3' scRNA-seq limitations

IsoCAPE is optimized for 10x Genomics 3' scRNA-seq (Chromium v2/v3). Limitations:

- **Sparse per-cell signal:** most genes have 1–5 UMI per cell; per-cell CE fraction is unreliable. Use pseudo-bulk (sum across cell type) for quantification.
- **3' bias:** only the terminal ~500bp of each transcript is captured; IsoCAPE cannot distinguish CE events near the annotated 3' end from normal APA.
- **Validation recommended:** bulk 3'-seq (QuantSeq, PAPERCLIP) or long-read sequencing (PacBio/Nanopore) should be used to confirm high-priority CE candidates.

---

## Repository Structure

```
IsoCAPE/
├── IsoCAPE/
│   └── scripts/
│       ├── bam_to_parquet.py           # Step 1: BAM → reads parquet
│       ├── bam_to_parquet_parallel.py  # Step 1: parallel multi-core version
│       ├── build_matrix.py             # Step 3: parquet → AnnData (.h5ad)
│       └── annotator/
│           ├── gtf_parser.py           # GTF/DB index builder
│           │   # intron_trees: protein_coding
│           │   # exon_trees: protein_coding + NMD + non_stop_decay
│           │   # known_3p_ends: protein_coding tx_biotype only
│           └── site_annotator.py       # Step 2: CE / PA / known annotation
│               # PAS_WINDOW=60bp, PAS_MIN_DIST=10bp
│               # known-window=50bp
│               # GN tag consistency filter
├── notebooks/
│   └── 01_AR-V7_demo.ipynb            # Validation: AR-V7 detection in CRPC
├── data/
│   └── cape_labels.csv                # Curated clinical label layer
├── results/
│   └── figures/
├── requirements.txt
└── README.md
```

---

## Results: Cryptic Polyadenylation in IDC Breast Cancer

IsoCAPE was applied to a public 10x Genomics 3' scRNA-seq dataset from invasive ductal carcinoma (IDC) breast cancer (5,680 cells, Cell Ranger v2, hg38). After QC and doublet removal, 2,974 cells were retained across 11 annotated cell types.

Candidate CE sites were first ranked by total read count and filtered by **tumor-versus-non-tumor** CE fraction differential (non-tumor reference: macrophages + stromal cells within the same dataset). Top candidates were further validated by IGV comparison against **healthy** PBMC and B cell/plasma cell BAM files. Two events satisfied all criteria for high-confidence cancer-specific cryptic polyadenylation:

---

### MXD4 — MYC antagonist inactivated by intronic polyadenylation

**MXD4** (MAX Dimerization Protein 4) is a transcriptional repressor and member of the MYC-MAX-MAD network. It antagonizes MYC oncogenic activity by competing for MAX binding, suppressing MYC-dependent transcriptional activation and cell transformation. Loss of MXD4 function — through mutation, deletion, or epigenetic silencing — relieves this brake on MYC, contributing to uncontrolled proliferation. MXD4 inactivation by intronic polyadenylation has been described in CLL (Lee et al., *Nature* 2018) and is now identified here in IDC breast cancer.

IsoCAPE detected two CE sites in **transcript intron 3** (chr4, − strand). Premature termination at this position retains only exons 1–3, producing a 3-exon isoform that lacks the bHLHZip domain required for MAX binding. This truncated protein cannot compete with MYC for MAX — functionally equivalent to loss of MXD4 — without any underlying DNA mutation.

---

### DNAAF1 — dynein assembly factor aberrantly expressed and truncated in tumor cells

**DNAAF1** (Dynein Axonemal Assembly Factor 1) orchestrates the cytoplasmic preassembly of dynein arm complexes required for ciliary motility. Its expression is normally restricted to ciliated epithelial cells and is absent from immune cells and most stromal populations. Loss-of-function mutations cause Primary Ciliary Dyskinesia (PCD), underscoring its non-redundant role in dynein biogenesis. In cancer, cilia loss is a hallmark of epithelial de-differentiation; the molecular drivers of this loss remain incompletely understood.

IsoCAPE detected four CE sites clustered in **transcript intron 9** (chr16, + strand), producing a truncated isoform that retains exons 1–9 but loses exons 10–13, including the C-terminal dynein assembly domain. Two features make this event particularly notable: DNAAF1 is aberrantly expressed in 61–67% of IDC tumor cells despite being absent from macrophages and stromal cells in the same dataset, and the CE event is confirmed tumor-specific by IGV. This constitutes a double hit — aberrant transcriptional activation paired with immediate isoform truncation — on a gene whose full-length product may normally constrain cell motility or vesicle trafficking.

---

### Figures

**Figure 1. MXD4 gene structure and IGV validation**

![MXD4 gene structure and IGV](results/figures/mxd4_structure_igv.png)

*Top: MXD4 gene structure (chr4, − strand, 6 exons). CE sites lie in transcript intron 3; the 3-exon truncated isoform loses the bHLHZip/MAX binding domain. IsoDecipher G1 (UTR=316bp, dominant) and G2 (UTR=1097bp) mark the annotated 3' ends. Bottom: IGV coverage tracks at the CE region. IDC breast cancer (red) shows a prominent pileup at both CE sites. PBMC (blue) and B cell/plasma cell (green) tracks show absent signal at the CE position despite MXD4 being expressed in these cell types — confirming tumor-specific intronic polyadenylation.*

---

**Figure 2. Single-cell UMAP: GEX vs CE signal**

![MXD4 DNAAF1 UMAP](results/figures/MDX4_DNAAF1_UMAP.png)

*MXD4 (top row): GEX is broadly distributed across all cell types, whereas CE signal is enriched in tumor cells and near-absent in macrophages and stromal cells — demonstrating that CE is a tumor-specific event independent of expression level. DNAAF1 (bottom row): GEX is aberrantly activated in tumor cells (absent in immune/stromal populations); CE signal mirrors this tumor-restricted pattern, concentrated in the same cells that aberrantly express the gene.*

---

**Figure 3. Pseudo-bulk CE signal by cell type**

![MXD4 DNAAF1 barplot](results/figures/MXD4_DNAAF1_barplot.png)

*Mean CE reads per cell across all annotated cell types (gray: non-tumor; red: tumor subtypes). MXD4 GEX is similar across cell types, yet MXD4 CE is 5–8× higher in tumor cells than in macrophages or stromal cells (Mann-Whitney U, p=1.75×10⁻¹⁰) — dissociating CE from expression level and implicating a tumor-specific splicing/polyadenylation mechanism. DNAAF1 CE is near-zero in all non-tumor cells and consistently elevated across every tumor subtype (p=2.21×10⁻¹³), consistent with the absence of DNAAF1 expression in immune and stromal populations.*

---

### Summary

| Gene | CE position | Intron | Non-tumor CE | Tumor CE | p-value |
|------|-------------|--------|-------------|----------|---------|
| MXD4 | chr4:2,254,615 (dominant) | transcript intron 3 | 0.009 | 0.141 | 1.75×10⁻¹⁰ |
| DNAAF1 | chr16:84,173,338 (dominant) | transcript intron 9 | 0.006 | 0.178 | 2.21×10⁻¹³ |

*CE signal = mean CE reads per cell. p-value: Mann-Whitney U, tumor vs non-tumor cells.*

Both events share a common theme: intronic polyadenylation as a mechanism to inactivate or truncate genes with potential tumor-suppressive or differentiation-associated functions, without requiring genetic mutation. MXD4 CE directly derepresses MYC — one of the most frequently activated oncogenes in breast cancer. DNAAF1 CE occurs on a background of aberrant gene activation, suggesting that the tumor transcriptional program both misexpresses and simultaneously truncates this gene. These findings illustrate the potential of IsoCAPE to detect functionally relevant isoform events from standard archived RNA-seq data, and motivate validation in larger cohorts and matched normal tissue.

---

- **Mouse / multi-species support** — mm10 GTF + reference genome
- **Spatial transcriptomics** — Visium BAM compatibility (spot-level APA resolution across tissue sections)
- **IsoCAPE-DB** — a community database of validated cryptic polyadenylation sites with clinical annotations across cancer types
- **IsoFormer integration** — IsoCAPE sites form the novel-site vocabulary for BAM-level foundation model training

---

## License

MIT License  
© 2026 Rene Yu-Hong Cheng
