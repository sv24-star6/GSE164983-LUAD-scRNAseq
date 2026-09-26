# Single-Cell Transcriptomic Characterisation of Human Lung Adenocarcinoma

A reproducible single-cell RNA-seq analysis of two primary human lung adenocarcinoma (LUAD) specimens from **GEO GSE164983**, implemented in R with Seurat.

The project covers quality control, sample-wise doublet detection, normalization, dimensionality reduction, unsupervised clustering, marker-guided cell-type annotation, tumour cellular-composition analysis, epithelial-associated gene analysis, and functional pathway enrichment.

## Project question

**What cellular populations and transcriptional programmes characterise the tumour microenvironment of primary human lung adenocarcinoma, and what biologically relevant pathways can be identified from single-cell RNA sequencing?**

## Dataset

- GEO accession: **GSE164983**
- Samples analysed: **2 primary human LUAD specimens**
- Technology: **10x Genomics Chromium Single Cell 3'**
- Input: processed filtered feature-barcode H5 matrices
- Reference genome: GRCh38

Raw sequencing files are not required for this workflow. Download the two processed `filtered_feature_bc_matrix.h5` files from GEO and place them in `data/`.

## Analysis workflow

```text
Processed 10x matrices
        |
        v
Initial dataset: 22,323 cells
        |
        v
Quality control
300 <= detected genes <= 6,000
mitochondrial RNA <= 20%
        |
        v
21,926 QC-passed cells
        |
        v
Sample-wise scDblFinder
2,517 predicted doublets removed
        |
        v
19,409 singlets
        |
        v
Log normalization + 2,000 HVGs
        |
        v
PCA (50 PCs)
        |
        v
Neighbour graph + clustering (30 PCs)
        |
        v
UMAP + marker-guided annotation
        |
        v
Cellular composition
        |
        v
Epithelial-associated marker analysis
        |
        v
GO / Reactome / KEGG enrichment
```

## Key results

### Quality control and doublet removal

The two samples contained **22,323 cells** before QC. Dataset-informed filtering retained **21,926 cells (98.2%)**. Doublets were detected independently within each 10x capture using scDblFinder. Predicted doublet rates were approximately **13.4%** in sample 591146 and **9.0%** in sample 614658, leaving **19,409 singlets** for downstream analysis.

### Cellular landscape

Marker-guided annotation identified broad immune, epithelial and stromal populations, including:

- CD4/conventional T cells
- Cytotoxic T cells
- Regulatory T cells
- Interferon-responsive T cells
- NK cells
- B cells and plasma cells
- Macrophage/monocyte cells
- Mast cells and neutrophils
- pDCs
- Epithelial cells
- Endothelial cells
- Fibroblasts

The two specimens showed substantial descriptive inter-tumour heterogeneity. Sample 591146 contained a more heterogeneous mixture with appreciable epithelial, cytotoxic T, NK and mast-cell populations, whereas sample 614658 was dominated by CD4/conventional T cells and contained a distinct interferon-responsive T-cell state.

### Epithelial-associated programme

Comparison of annotated epithelial cells with the other cell populations identified **1,880 epithelial-associated marker genes**. Strongly enriched genes included `EHF`, `CXCL17`, `TFF3`, `CEACAM6`, `PIGR`, `ITGB6`, `SLC6A14` and `BPIFB1`.

Functional enrichment showed convergent biological themes involving:

- cell-junction organisation
- extracellular-matrix organisation and ECM-receptor interactions
- integrin and focal-adhesion signalling
- tight and adherens junctions
- RHO-family GTPase signalling
- receptor tyrosine kinase signalling

These results describe programmes enriched in epithelial cells relative to other cell populations in the dataset; they are **not** tumour-versus-normal differential-expression results.

## Why the samples were not forcibly integrated

UMAP and PCA showed sample-associated structure. With only two patients, biological inter-patient variation and technical sample effects cannot be reliably disentangled. Therefore, the primary analysis preserves the unintegrated structure rather than automatically applying Harmony/CCA and potentially removing genuine biological variation.

## Repository structure

```text
GSE164983-LUAD-scRNAseq/
├── R/
│   ├── 01_load_data.R
│   ├── 02_quality_control.R
│   ├── 03_doublet_detection.R
│   ├── 04_normalization_clustering.R
│   ├── 05_cell_type_annotation.R
│   ├── 06_epithelial_analysis.R
│   └── 07_pathway_enrichment.R
├── data/       # GEO H5 files; excluded from Git
├── figures/    # Generated analysis figures
├── results/    # Lightweight result tables
├── .gitignore
└── README.md
```

## Reproducing the analysis

Run the scripts sequentially from the repository root:

```text
01_load_data.R
02_quality_control.R
03_doublet_detection.R
04_normalization_clustering.R
05_cell_type_annotation.R
06_epithelial_analysis.R
07_pathway_enrichment.R
```

Main R packages include **Seurat**, **SingleCellExperiment**, **scDblFinder**, **gprofiler2**, **dplyr**, **ggplot2**, **patchwork**, and **hdf5r**.

## Interpretation and limitations

This is an exploratory analysis of **two LUAD specimens**, not a population-level cohort study. Sample identity is confounded with patient-specific biology and potential technical effects, so differences between the two samples are presented descriptively rather than as generalisable LUAD effects. The dataset does not contain matched normal controls. Epithelial cells are therefore labelled **epithelial**, not automatically classified as malignant; establishing malignancy would require additional evidence such as copy-number inference or appropriate reference cells.

The epithelial enrichment analysis compares epithelial cells with other annotated cell populations. Consequently, enriched pathways should not be interpreted as pathways upregulated in LUAD relative to normal lung.

## Skills demonstrated

Single-cell RNA sequencing · Seurat · scDblFinder · quality control · dimensionality reduction · unsupervised clustering · marker-based cell annotation · differential marker analysis · functional enrichment · GO · Reactome · KEGG · g:Profiler · reproducible R workflows · Git/GitHub

## Data availability

The source data are publicly available through **NCBI Gene Expression Omnibus, accession GSE164983**. Large raw/processed data files and intermediate Seurat `.rds` objects are intentionally excluded from this repository.
