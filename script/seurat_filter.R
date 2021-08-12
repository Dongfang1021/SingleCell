
# R script used for cell filteration
.libPaths("/ifs/TJPROJ3/HW/PUBLIC/software/Rlib/Rlib")
library(argparse)
#library(dplyr)
library(Seurat)
#library(reshape2)
library(ggplot2)

parser = ArgumentParser(description="Seurat analysis for sc-RNAseq")
parser$add_argument("--input", help="input file or directory. Required.", required=T)
parser$add_argument("--input_format", help="input format, default: h5", choices=c("h5","csvd","10x","csv",'rda'), default="csvd")
parser$add_argument("--outdir", help="output directory. Required.", required=T)
parser$add_argument("--prefix", help="sample name or group name. Required.", required=T)
qc_group = parser$add_argument_group("QC group", "arguments for QC")
qc_group$add_argument("--min_cells", help="Filter genes, genes expressed in min cells", type='integer', default=3)
qc_group$add_argument("--low_nGene", help="cells have low genes", type='integer', default=200)
qc_group$add_argument("--high_nGene", help="cells have high genes", type='integer', default=8000)
qc_group$add_argument("--high_pMT", help="cells have high percent.mito",type='integer', default=50)
qc_group$add_argument("--rRNA_genes", help="rRNA genelist, default: NULL")
qc_group$add_argument("--genome",help="genome version", choices=c("GRCh38","hg19","mm10"),required=T)
args = parser$parse_args()
str(args)

input = args$input
input_format = args$input_format
outdir = args$outdir
prefix = args$prefix
genome = args$genome
if(genome == "GRCh38" || genome == "hg19"){
	pattern_MT="^MT-"
}
if(genome == "mm10"){
	pattern_MT="^mt-"
}
## qc
min_cells = args$min_cells
low_feature = args$low_nGene
high_feature = args$high_nGene
high_mt_percent = args$high_pMT

data_10X <- Read10X(data.dir = input, gene.column = 2, unique.features = TRUE)


#nUMI=nCount_RNA；nGene=nFeature_RNA
barcode_matrix_raw <- CreateSeuratObject(counts = data_10X,project = "10X")
barcode_matrix_raw[["percent.MT"]] <- PercentageFeatureSet(barcode_matrix_raw, pattern = pattern_MT)
barcode_matrix_raw1<-cbind(rownames(barcode_matrix_raw@meta.data),barcode_matrix_raw@meta.data)
colnames(barcode_matrix_raw1)[1]<-"Barcode"
write.table(barcode_matrix_raw1[c(1,3,4,5)], file.path(outdir, paste(prefix,'_BarcodeMatrix_meta.xls',sep='')),sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)

#plot function dual.plot
low_res = 70
mid_res = 150
high_res = 300
dual.plot <- function(fig, file.prefix, w=7, h=7, res=75){
        pdf(paste(file.prefix,".pdf",sep=""), width = w, height = h)
        print(fig)
        dev.off()
        png(paste(file.prefix,".png",sep=""), width = w*res, height = h*res, res = res, type="cairo-png")
        print(fig)
        dev.off()
}
#plot
feature.names = c("nFeature_RNA","nCount_RNA","percent.MT")
scatter_plots = list()
for (i in 2:length(feature.names)){
        scatter_plots[[i-1]] = FeatureScatter(object = barcode_matrix_raw, feature1 = feature.names[1], feature2 = feature.names[i], group.by="orig.ident")
}
nscatter = length(scatter_plots)
f = CombinePlots(plots = scatter_plots, ncol = nscatter, legend = "none")
dual.plot(f, file.path(outdir, paste(prefix,'_BarcodeMatrix_meta.scatter',sep="")), w=nscatter*5, h=5,res=mid_res)

p = VlnPlot(barcode_matrix_raw, features = feature.names, ncol = length(feature.names), group.by="orig.ident")
dual.plot(p, file.path(outdir, paste(prefix,'_BarcodeMatrix_meta.violin',sep="")), w=4*length(feature.names),h=6,res=mid_res)

#filter
barcode_matrix <- CreateSeuratObject(counts = data_10X,min.cells = min_cells,min.features = low_feature,project = prefix)
barcode_matrix[["percent.MT"]] <- PercentageFeatureSet(barcode_matrix, pattern = pattern_MT)
#barcode_matrix <- WhichCells(barcode_matrix, expression = percent.MT < 50)
barcode_matrix <- subset(x =barcode_matrix, subset=nFeature_RNA < high_feature & percent.MT < high_mt_percent)
barcode_matrix2<-cbind(rownames(barcode_matrix@meta.data),barcode_matrix@meta.data)
colnames(barcode_matrix2)[1]<-"Barcode"
write.table(barcode_matrix2[c(1,3,4,5)], file.path(outdir, paste(prefix,"_BarcodeMatrix_filter.xls",sep='')),sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)

#plot 
scatter_plots = list()
for (i in 2:length(feature.names)){
        scatter_plots[[i-1]] = FeatureScatter(object = barcode_matrix, feature1 = feature.names[1], feature2 = feature.names[i], group.by="orig.ident")
}
nscatter = length(scatter_plots)
f = CombinePlots(plots = scatter_plots, ncol = nscatter, legend = "none")
dual.plot(f, file.path(outdir, paste(prefix,'_BarcodeMatrix_filter.scatter',sep="")), w=nscatter*5, h=5,res=mid_res)

p = VlnPlot(barcode_matrix, features = feature.names, ncol = length(feature.names), group.by="orig.ident")
dual.plot(p, file.path(outdir, paste(prefix,'_BarcodeMatrix_filter.violin',sep="")), w=4*length(feature.names),h=6,res=mid_res)
