# Session Data Overview

## Beginner Sessions
- **Monday - Introduction to single-cell and Seurat – through making a Seurat object**
    - Input:
        - `input/scRNA/GSM6925133_OX1X/` files
          (e.g. the three matrix files from the GSM6925133_OX1X dataset)
        - Optionally an h5 file that bundles the matrices
    - Output:
        - `output/rdata/OX1X.qs`
          (a Seurat object containing the single-cell data from the GSM6925133_OX1X dataset)

- **Tuesday - Processing single-cell data and cell type annotation**
    - Input:
        - `output/rdata_premade/OX1X.qs` (or the attendee’s object from Monday)
    - Output:
        - `output/rdata/OX1X_annotated.qs`
          (a Seurat object from the GSM6925133_OX1X dataset with cell type annotations and processing applied)

- **Wednesday - Combining datasets and data integration**
    - Input:
        - `output/rdata/OX1X_annotated.qs` (or the attendee’s object from Tuesday)
        - `input/rdata_premade/scRNA_other_annotated.qs`
          (a premade list of annotated brain Seurat objects from additional datasets)
    - Output:
        - `output/rdata/brain_scRNA_combined.qs`
          (a single Seurat object resulting from combining the GSM6925133_OX1X dataset with the other brain datasets)

- **Thursday - Differential expression analysis, GSEA, and experimental design**
    - Input:
        - `output/rdata_premade/brain_scRNA_combined.qs` (or the attendee’s combined object from Wednesday)
    - Output:
        - Processed results (tables and plots) saved in `output/rdata/`
          (including differential expression CSV files, volcano plots, and annotated Seurat objects)

## Advanced Sessions
- **Monday - Understanding Seurat Object Structure**
    - Input:
        - `output/rdata/OX1X_annotated.qs`
          (a premade Seurat object designed for “under the hood” exploration)
    - Output:
        - none

- **Tuesday - Spatial Single-Cell Analysis**
    - Input:
        - `input/spatial/` data
          (a Visium spatial transcriptomics dataset from a public resource)
    - Output:
        - none

- **Wednesday - Cell-Cell Communication Analysis**
    - Input:
        - `output/rdata/brain_scRNA_combined.qs`
          (the integrated brain Seurat object from the beginner sessions, further processed for cell–cell communication)
    - Output:
        - none

- **Thursday - 10X Multiomics Analysis**
    - Input:
        - `input/multiomics/` data
          (a raw multiomic dataset from 10X Multiome runs, including raw GEX and ATAC matrices)
    - Output:
        - none