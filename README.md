# Seurat scRNA-seq Analysis, Powered by BadranSeq & CyteTypeR

**The Seurat PBMC3k tutorial, rebuilt with publication-ready figures and automated cell-type annotation.**

- Every visualization rendered **twice** — Seurat default, then the [BadranSeq](https://github.com/wolf5996/BadranSeq) alternative — so you can see the difference yourself
- Cell-type annotation done **twice** — manually with canonical markers, then automatically with [CyteTypeR](https://github.com/NygenAnalytics/CyteTypeR) — zero effort, full comparison
- Seurat's defaults get the job done. BadranSeq makes them belong in a paper.

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

---

## The Good Stuff

**Annotated UMAP** — cell borders, white-fill labels, vivid palette. `BadranSeq::do_UmapPlot()` doing what `DimPlot()` wishes it could:

<img src="readme_figures/annotated_umap_badranseq.png" alt="Annotated UMAP with BadranSeq" width="100%" />

**Marker gene expression** — nine canonical markers, viridis scaling, cell borders. `BadranSeq::do_FeaturePlot()`:

<img src="readme_figures/marker_features_badranseq.png" alt="Marker feature plots with BadranSeq" width="100%" />

**Manual vs CyteTypeR annotation** — hand-curated labels vs API-driven predictions, side by side:

*Manual annotation (canonical markers + prior knowledge):*

<img src="readme_figures/annotation_manual.png" alt="Manual cell type annotation" width="100%" />

*CyteTypeR automated annotation:*

<img src="readme_figures/annotation_cytetype.png" alt="CyteTypeR automated annotation" width="100%" />

**Statistical violin plots** — Kruskal-Wallis omnibus tests, median annotations, jittered points. No more naked violins:

<img src="readme_figures/stats_violin_badranseq.png" alt="Statistical violin plots with BadranSeq" width="100%" />

---

## Why BadranSeq?

| Seurat Default | BadranSeq Replacement | What You Get |
|---|---|---|
| `DimPlot()` | `do_UmapPlot()` | Cell borders, auto cluster labels, vivid palette |
| `DimPlot(reduction = "pca")` | `do_PcaPlot()` | Variance-explained % on axes |
| `FeaturePlot()` | `do_FeaturePlot()` | Viridis scaling, cell borders, cleaner layout |
| `VlnPlot()` | `do_ViolinPlot()` | Boxplot overlay, median line, jittered points |
| `ElbowPlot()` | `EnhancedElbowPlot()` | Variance explained %, cutoff line |
| N/A | `fetch_cell_data()` | Tidy metadata for custom ggplot2 QC plots |
| Manual `RenameIdents()` | **CyteTypeR** | Automated annotation + interactive report |

---

## Quick Start

```r
# Dependencies
install.packages(c("Seurat", "remotes"))
remotes::install_github("satijalab/seurat-data")
remotes::install_github("wolf5996/BadranSeq")
remotes::install_github("NygenAnalytics/CyteTypeR")
```

```bash
# Render
quarto render seurat_analysis_powered_by_badranseq_guide.qmd
```

---

## Project Layout

```
seurat_analysis_powered_by_badranseq/
├── docs/              # Plans and documentation
├── read/              # Input data
├── scripts/           # Analysis code (git root)
│   ├── seurat_analysis_powered_by_badranseq_guide.qmd
│   ├── readme_figures/ # README showcase images
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
| CyteTypeR Interactive Report | [View Report](https://nygen-labs-prod--cytetype-api.modal.run/report/34fac9e9-3c43-4c46-95f4-6b2994e57ada) |

---

**Author:** Badran Elshenawy — University of Oxford

Made with coffee and an unreasonable number of UMAP plots.
