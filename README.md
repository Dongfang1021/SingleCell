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
Cell Ranger is a set of analysis pipelines that process Chromium single-cell data to align reads, generate feature-barcode matrices, perform clustering and other secondary analysis, and more. Cell Ranger includes four pipelines relevant to the 3' Single Cell Gene Expression Solution and related products.
### 2.1 From BCL to Fastq
`cellranger mkfastq` is used to demulitiplex BCL files (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_fq). It is a wrapper around Illumina's bcl2fastq, with additional features that are specific to 10X libraries and a simplified sample sheet format.
```shell
cellranger mkfastq --id={output_foldername} \
		   --run={flowcell_runid_path} \
		   --csv=samplesheet.csv
```
Final fragment from Chromium single cell 3' v3 library
![](image/lib-v3.png)
- I1.fastq.gz: 8bp P7 index (sample index)
- R1.fastq.gz: total 28 bp, 16bp barcode + 10bp UMI
- R2.fastq.gz: 91bp cDNA
  
The FASTQ files are named according to the sample column of the sample sheet. If a sample ID was not specified, the flowcell ID is used instead(not shown here). In addition to the FASTQ files, bcl2fastq generates various summary files. If -stats-dir was not specified, summary and statistic files will be stored in a Stats folder by default.


### 2.2 Data Process and Quality Control
Reads from scRNA-seq were processed by Cell Ranger software (10X Genomics) with the default parameters for each sample seperately. In short, alignment to the reference, UMI(unique molecular identifier) collapsing. UMI counting and initial quality control was performed.The filtered gene expression matrices containing only celluar barcodes were generated in this step. Then Seurat R package was further used to get data quality control.
#### 2.2.1 Cellranger count
`Cellranger count` quantifies single-cell gene expression. Meanwhile, it provides sequencing summary
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct
```shell
cellranger count --id={sample} --transcriptome={reference genome} --fastqs={fastq files} 
```
##### Initial quality control
```python
python qc.py --input --output
```
##### Alignment Summary
```python
python alignment.py --input --output
```
##### Cell Calling and UMI Counting
```python
python cell_umi.py --input --output
```


#### 2.2.2 Quality control and cell filtering based on Seurat R package
The "filtered_feature_bc_matrix" generated by Cell Ranger were read into the Seurat R package. For each cell, quality control metrics such as the total number of counts and the proportion of counts in mitochondrial genes were calculated.
```R
Rscript qc_seurat.R --input --output
```
Cells that met any one of the following criteria were filtered out for downstream processing in each sample: < 200 genes with present in each cell, genes with non-zero counts in at most 3 cells, > 50,000 feature count(potential multiplets), or the proportion of the feature count attributable to mitochondrial genes was greater than 50%
seurat is used for filtering
```R
Rscript filter.R --input --output
```

### 2.3 Identification of Highly Variable Genes
Feature selection removes the uniformative genes and identifies the most relevant features to reduce the number of dimensions used in downstream analysis. The identification of highly variable genes (HVGs) was based on the assumption that the highly variable expression of these genes across cells is resluted from biological effects rather than merely technical noise.
After cell filtering, 10X gene expression matrices from each sample were loaded into Seurat. The expressions of each gene were normalized by total counts for that cell, mulitplied by a scale factor (the median UMI counts for all cells within the sample), and natural-log transformed using log1p. Then the expressions of each gene were scaled.
Next, highly variable genes were identified using FindVariableGenes function with the highest standardized variance selected by selection method = 'vst'. The top 3000 most variable genes of each sample selected by Seurat were used to compute the PCs (principal components).
Normalized expression of cells
```R
Rscript normalization.R --input --output
```

Highly variable genes were identified using FindVariableGenes function
```R
Rscript HVG.R --input --output
```
### 2.4 Cell subpopulation identification
A key goal of 10X scRNA-seq data analysis is to identify cell subpopulations (different populations are often distinct cell types) within a certain condition or tissue to unravel the heterogeneity of cells. To identify a gene expression signature associated with this sample or group of cells, HVGs previously determined are used as input for dimensionality reduction via principal component analysis (PCA). The resulting PCs were then used as input for clustering analysis.
#### 2.4.1 Principal Components Analysis

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
