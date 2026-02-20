# Seurat scRNA-seq Analysis, Powered by BadranSeq

**A complete PBMC3k single-cell RNA-seq walkthrough — every plot shown twice: Seurat default, then the BadranSeq alternative.**

> **TL;DR** — Follow the canonical Seurat PBMC3k tutorial from raw counts to annotated clusters.
> At every visualisation step, compare Seurat's built-in plots against BadranSeq replacements
> (better palettes, cell borders, axis labels that actually mean something).
> Cell-type annotation is handled automatically by CyteTypeR, then compared head-to-head with manual labels.

---

## Pipeline Overview

```
PBMC3k Data                                                              Annotate
    │                                                                       ▲
    ▼                                                                       │
   QC ──→ Normalize ──→ Variable Features ──→ Scale ──→ PCA ──→ Cluster ──→ UMAP ──→ DEA ──┤
    │                                                     │        │         │               │
    │                                                     │        │         │               │
 BadranSeq                                             BadranSeq  BadranSeq BadranSeq    CyteTypeR
 fetch_cell_data()                                     do_PcaPlot  Enhanced  do_UmapPlot  Automated
 QC histograms                                         ElbowPlot   DimPlot  do_FeaturePlot Annotation
                                                                            do_ViolinPlot
```

Every step that produces a figure renders it **twice** — once with the Seurat default, once with the BadranSeq replacement — so you can compare them side by side.

---

## Why BadranSeq?

| Seurat Default | BadranSeq Replacement | Why Bother |
|---|---|---|
| `DimPlot()` | `do_UmapPlot()` | Cell borders, automatic cluster labels, vivid colour palette |
| `DimPlot(reduction = "pca")` | `do_PcaPlot()` | Auto variance-explained percentages on axes |
| `FeaturePlot()` | `do_FeaturePlot()` | Viridis colour scaling, cell borders, cleaner layout |
| `VlnPlot()` | `do_ViolinPlot()` | Boxplot overlay, median line, jittered points |
| `ElbowPlot()` | `EnhancedElbowPlot()` | Cumulative variance explained, cutoff indicator line |
| N/A | `fetch_cell_data()` | Tidy metadata extraction for custom ggplot2 QC histograms |
| Manual `RenameIdents()` | **CyteTypeR** | API-driven automated cell-type annotation + interactive HTML report |

---

## Quick Start

### 1. Install dependencies

```r
# Install Seurat and SeuratData
install.packages("Seurat")
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")

# Install BadranSeq
remotes::install_github("wolf5996/BadranSeq")

# Install CyteTypeR
remotes::install_github("NygenAnalytics/CyteTypeR")
```

### 2. Render the guide

```bash
quarto render seurat_analysis_powered_by_badranseq_guide.qmd
```

The rendered HTML will appear in the same directory.

---

## Project Layout

```
seurat_analysis_powered_by_badranseq/
├── docs/              # Plans and documentation
├── read/              # Input data
├── scripts/           # Analysis code (git root)
│   ├── seurat_analysis_powered_by_badranseq_guide.qmd
│   └── scripts.Rproj
├── checkpoints/       # Intermediate .rds files
└── write/
    ├── figures/        # All output figures
    └── tables/         # Marker gene tables
```

---

## Requirements

| Dependency | Version |
|---|---|
| R | >= 4.0 |
| Quarto | >= 1.3 |
| Seurat | >= 5.0 |
| SeuratData | latest |
| BadranSeq | latest (GitHub) |
| CyteTypeR | latest (GitHub) |

---

## Resources

| Resource | Link |
|---|---|
| BadranSeq | <https://github.com/wolf5996/BadranSeq> |
| CyteTypeR | <https://github.com/NygenAnalytics/CyteTypeR> |
| Seurat PBMC3k Tutorial | <https://satijalab.org/seurat/articles/pbmc3k_tutorial> |
| CyteTypeR Interactive Report | [View Report](https://nygen-labs-prod--cytetype-api.modal.run/report/35259c0e-5d4d-47a1-a128-69861ba5d08c) |

---

**Author:** Badran Elshenawy — University of Oxford

Made with coffee and an unreasonable number of UMAP plots.
