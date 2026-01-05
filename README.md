# AlbaoRunx3Manuscript

Code and analysis assets for the manuscript investigating RUNX3-dependent regulation in CD8 T cells. The repository collects single-cell RNA-seq processing notebooks, ATAC-seq differential accessibility workflows, figure exports, and helper scripts used throughout the study.

## Repository layout

- `single_cell/` — Scanpy-based notebooks and outputs for preprocessing, clustering, differential expression, pathway enrichment, signature scoring, and RNA velocity analyses.
- `atac/` — ATAC-seq peak filtering, differential accessibility, peak set enrichment analysis (PSEA), and motif analyses implemented in Jupyter notebooks and exported HTML reports.
- `csv/` — Compressed intermediate tables for clustering, differential expression, GSEA outputs, and figure-ready summaries.
- `figures/` — Generated figure assets grouped by panel.
- `reanalysis/` — Supplemental GSEA notebook for comparative datasets.
- `scripts/` — Utility R scripts for flow cytometry processing, SEA plotting, and notebook conversion.

## Analysis overview

### Single-cell RNA sequencing

- P14 CD8 T cells were hash tagged with BioLegend TotalSeq A antibodies (A0301–A0309) before pooling and loading 100,000 multiplexed cells across two 10x Genomics Single Cell 3' v3.1 GEMs.
- Hash tag oligo (HTO) libraries were amplified from cDNA with the BioLegend HTO additive primer, while gene expression (GEX) libraries followed the standard 10x protocol. Sequencing depth targeted ~35,000 reads per cell for GEX and ~300 reads per cell for HTO on an Illumina NextSeq 2000.
- CellRanger v7.1.0 with the 2020-A mm10 reference generated count matrices. Quality control in Scanpy v1.9.5 included SoupX ambient RNA correction, Scrublet and scDblFinder doublet removal, and filtering on mitochondrial content (>5%), UMI count (<3,000), and detected genes (<1,250). Cell cycle effects were regressed using S and G2/M gene lists, and Leiden clustering at resolution 1.0 (14–16 clusters) guided downstream analyses.
- Differential expression relied on the Scanpy Wilcoxon rank-sum test (v1.9.6), with GSEA run on Wilcoxon statistics via clusterProfiler v4.14.0 and DOSE v4.0.0 (R v4.4.2). Transcriptome correlations used Spearman coefficients on mean-normalized counts of the top 2,000 variable genes. RNA velocity inputs were produced with velocyto v0.17 and analyzed with scVelo v0.3.1.

### ATAC-seq

- Nuclei from sorted P14 CD8 T cells were lysed and permeabilized for in situ Tn5 transposition (Nextera DNA Library Preparation kit), titrated to 1.25 μL enzyme per 5×10⁴ nuclei in a 50 μL reaction. Libraries were PCR-amplified to optimal cycles, quality-checked on an Agilent TapeStation (targeting a 2:1 170–280 bp to 280–390 bp fragment ratio), and sequenced on an Illumina NextSeq 2000 at ~0.9× whole-genome coverage (paired-end 50×2).
- Processing and peak calling used the nf-core/atacseq v2.1.2 pipeline (Nextflow 24.04.2) against mm10 in broad-peak mode. Public datasets (e.g., GSE111149, GSE88987, GSE144383, GSE131871, GSE213041) were merged with study data to define consensus peaks; peaks present in at least one replicate from three datasets yielded 47,192 shared peaks from an initial 285,385. Differential accessibility employed DESeq2 v1.46.0, PSEA used clusterProfiler/DOSE, and motif analysis relied on HOMER v5.1. Peak overlaps with ChIP-seq datasets (RUNX3, TBET, c-JUN) were computed with bedtools v2.31.1 (`intersect -f 0.5 -F 0.5`), with accessibility profiles generated via deepTools v3.5.6.

### Gene and peak resources

- CD8 T cell-specific gene lists follow prior definitions, supplemented by re-analyses of public bulk RNA-seq datasets (e.g., Tcf7-reporter, ex-KLRG1, Tox genotypes, terminal-TEM). Reads were mapped with Salmon v1.10.2 (mm39, Ensembl 110) using dataset-specific k-mer parameters, with differential expression performed via DESeq2 v1.46.0 on TPM values.
- ATAC-seq peak lists were generated ad hoc from the consensus procedure above to support PSEA.

### Computational reproducibility

Containerized tools were used end to end:

- Public images:
  - CellRanger v7.1.0 (`docker.io/litd/docker-cellranger:v7.1.0`) for alignment and count matrix generation.
  - bedtools v2.31.1 (`docker.io/staphb/bedtools:2.31.1`) for peak overlaps.
  - deepTools v3.5.6 (`quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0`) for accessibility profile generation.
  - Salmon v1.10.2 (`docker.io/combinelab/salmon:1.10.2`) for transcript quantification.
- Custom DockerHub images:
  - HOMER v5.1 (`docker.io/pipkinlab/homer:5.1`) for motif analysis.
  - Scanpy v1.9.5 (`docker.io/pipkinlab/scanpy:1.9.5`) for core scRNA-seq QC and clustering.
  - Scanpy v1.9.6 with scVelo (`docker.io/pipkinlab/scanpy:1.9.6`) for differential expression and RNA velocity.
  - R 4.4.2 with DESeq2 v1.46.0 and plotting stack (`docker.io/pipkinlab/r-deseq2-plotting:4.4.2`) for ATAC-seq differential accessibility and visualization.

Gene and peak sets are available at [ScrippsPipkinLab/CommonGeneSets](https://github.com/ScrippsPipkinLab/CommonGeneSets). Code for this study is hosted in this repository.
