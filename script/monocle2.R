###
# Author Dongfang Hu
# This script is designed to pipeline Seurat object output into Monocle2
# It will start from reading in a Seurat object with the count sparse matrix, UMAP coordinates, and cluster information
# Date 20200701

####

library(Seurat)
library(monocle)
library(argparse)
parser = ArgumentParser(description="Seurat analysis for sc-RNAseq")
parser$add_argument("--input", help="input files or directories. Required.", required=T)
parser$add_argument("--outdir", help="analysis directory. Required.", required=T)
parser$add_argument("--prefix", help="sample or group",required=T)
args = parser$parse_args()
str(args)

inputfile = args$input
outdir = args$outdir
prefix = args$prefix

dir.create(outdir,showWarnings = FALSE,recursive = TRUE)

#function
dual.plot2 <- function(fig, file.prefix, w=7, h=7, res=75){
        pdf(paste(file.prefix,".pdf",sep=""), width = w, height = h)
        print(fig)
        dev.off()
        png(paste(file.prefix,".png",sep=""), width = w*res, height = h*res, res = res, type="cairo-png")
        print(fig)
        dev.off()
}

### Reading in Seurat object
print("Readingin Seurat objects...")
seurat <- readRDS(file = inputfile)

#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(seurat@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
print("Construct monocle newCellDataSet...")
cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              #lowerDetectionLimit = 0.5,
                              expressionFamily = uninormal())# since I have already normalized, thresholded and scalled in Suerat v3.

#Run ordering algorithm
var_genes <- seurat[["RNA"]]@var.features
ordering_genes <- var_genes
cds <- setOrderingFilter(cds, ordering_genes)
print(dim(exprs(cds)))

## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
cds <- reduceDimension(cds,norm_method="none", 
                        reduction_method="DDRTree",
                        max_components=4,
                        scaling=TRUE,
                        verbose=TRUE,
                        pseudo_expr=0)

# First decide what you want to color your cells by
print(head(pData(cds)))

## order cells change colors and theta to match your plot
cds <- orderCells(cds)
write.table(pData(cds), file = file.path(outdir,paste(prefix,".Trajectory_BarcodeMatrix.xls",sep="")), sep='\t', quote=F, row.names=T)

#plot
p <- plot_cell_trajectory(cds, 
                     color_by = "seurat_clusters",
                     theta = -15,
                     show_tree = TRUE)
dual.plot2(p, file.path(outdir,paste(prefix,".TrajectoryPlot_by_clusters",sep="")),w=8, h=5,res=150)

p <- plot_cell_trajectory(cds, color_by = "Pseudotime",theta = -15,show_tree = TRUE)
dual.plot2(p, file.path(outdir,paste(prefix,".TrajectoryPlot_by_Pseudotime",sep="")),w=8, h=5,res=150)




my_genes <- head(row.names(subset(fData(cds))))
cds_subset <- cds[my_genes, ]
f <- plot_genes_in_pseudotime(cds_subset, color_by = "orig.ident")
dual.plot2(f, file.path(outdir,paste(prefix,"Plot_genes_in_Pseudotime",sep="")),w=8, h=5,res=150)


