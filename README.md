# 10X Single Cell RNAseq workflow
## 1. Introduction to single-cell RNA-seq
### 1.1 Bulk RNA-seq
Bulk RNAseq technologies have been widely used to study gene expression patterns at population level in the past decade.
- A major breakthrough(replacement of microarray) in the late 00's and has been widely used since.
- Measures the average expression level for each gene across a large population of input cells.
- Useful for comparative transcriptomics, e.g. samples of the same tissue from different species
- Useful for quantifying expression signatures from ensembles, e.g. in disease studies
- Insufficient for studying heterogeneous systems, e.g. early development studies, complex tissues(brain)
- Does not provide insights into the stochastic nature of gene expression
### 1.2 scRNA-seq
 The advent of single-cell RNA sequencing (scRNA-seq) provides unprecedented opportunities for exploring gene expression profile at the single-cell level. Currently, scRNA-seq has become a favorable choice for studying the key biological questions of cell heterogeneity and the development of early embryos (only include a few number of cells), since bulk RNA-seq mainly reflects the averaged gene expression across thousands of cells.The recently developed droplet-based single-cell transcriptome sequencing (scRNA-seq) technology makes it feasible to perform a population-scale scRNA-seq study, in which the transcriptome is measured for tenes of thousands of single cells from multiple individuals.
 - A new technology, first publiction by Tang et al, 2009.
 - Did not gain widespread popularity until 2014 when new protocols and lower sequencing costs made it more accessible
 - Measures the distribution of expression level for each gene across a population of cells
 - Allows to study new biological questions in which cell-specific changes in transcriptome are important, e.g. cell type identification, heterogeneity of cell response, stochasticity of gene expression, inference of gene regulatory networks across the cells.
 - Datasets range from 100 to 1000000 cells and increase size every year.
 - Currently there are several different protocols in use, e.g. SMART-seq2, CELL-seq and Drop-seq.
  ![](/image/single_cell_lib.jpeg)
  Figure 1.1 different protocols to prep single cell library
  Hwang, B., Lee, J.H. & Bang, D. Single-cell RNA sequencing technologies and bioinformatics pipelines. Exp Mol Med 50, 1–14 (2018).
 - There are also commercial platforms available, including the Fluidigm C1, Wafergen ICELL8 and the 10X Genomics Chromium
 - Several computational analysis methods from bulk RNA-seq can be used
 - In most cases computational analysis requires adaption of the existing methods or development of new ones.
Note: got lots of idea from website "Analysis of single cell RNA-seq data" https://scrnaseq-course.cog.sanger.ac.uk/website/index.html
## 2. Data analysis Workflow
![](/image/10X_pipeline.png)
### 2.1 From BCL to Fastq
`cellranger mkfastq` is used to demulitiplex BCL files (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_fq)
```shell
cellranger mkfastq --id={output_foldername} \
		   --run={flowcell_runid_path} \
		   --csv=samplesheet.csv
```
![](image/lib-v3.png)
Data Processing and Quality Control

### 2.2 Count summary
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct
```shell
cellranger count --id={sample} --transcriptome={reference genome} --fastqs={fastq files} 
```
#### Alignment Summary
Cellranger uses an aligner called STAR, which performs splicing-aware alignment of reads to genome.

#### Cell Calling and UMI Counting
#### Quality control and Cell filtering
seurate is used for quality control
```R
```
seurat is used for filtering
```R
```

### 2.3 Identification of Highly Variable Genes
seurat: FindVariable Genes
```R
```
### 2.4 Cell subpopulation identification
seurat: FindClusters
```R
```
### 2.5 Marker Gene Detection
seurat: FindAllMarkers
```R
```
### 2.6 Enrichment analysis
GO, KEGG, REACTOME, GSEA
`clusterprofile` is used for enrichment analysis 
```R
```
### 2.7 Functional annotation of transcription factor
`TFCat` is used for functional annotation of transcription factor.
### 2.8 Protein-protein interaction network analysis
STRING
