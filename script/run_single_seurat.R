# Rscript used for scRNA analysis
# This script covers most of scRNA analysis

library(argparse)
library(dplyr)
library(Seurat)
library(umap)
library(reshape2)
library(ggplot2)
library(psych)
library(pheatmap)

parser = ArgumentParser(description="Seurat analysis for sc-RNAseq")
parser$add_argument("--input", help="input file or directory. Required.", required=T)
parser$add_argument("--input-format", help="input format, default: csvd", choices=c("h5","csvd","10x","csv",'rda'), default="csvd")
parser$add_argument("--prefix", help="prefix name or group name. Required.", required=T)
parser$add_argument("--gene-name", help="gene name annotation file")
parser$add_argument("--analydir", help="analysis directory. Required.", required=T)
qc_group = parser$add_argument_group("QC group", "arguments for QC")
qc_group$add_argument("--min-cells", help="Filter genes, genes expressed in min cells", type='integer', default=3)
qc_group$add_argument("--low-nGene", help="cells have low genes", type='integer', default=200)
qc_group$add_argument("--high-nGene", help="cells have high genes", type='integer', default=8000)
qc_group$add_argument("--high-pMT", help="cells have high percent.mito",type='integer', default=50)
qc_group$add_argument("--genome",help="genome version", choices=c("GRCh38","hg19","mm10"),required=T)
qc_group$add_argument("--rRNA-genes", help="rRNA genelist, default: NULL")
qc_group$add_argument("--normalize-method", help="normalize method, default: LogNormalize", choices=c('LogNormalize','CLR','RC','SCT'), default='LogNormalize')
pca_group = parser$add_argument_group("PCA group", "arguments for PCA analysis")
pca_group$add_argument("--selection-method", help="selection method for FindVariableFeatures, default: logNorm", choices=c('mvp','vst','disp'), default='vst')
pca_group$add_argument("--dispersion-cutoff", help="low cutoff for feature dispersions, default: 1", type='double', default=1)
pca_group$add_argument("--feature-number", help="number of VariableFeatures for PCA, default: 3000", type='integer', default=3000)
pca_group$add_argument("--pc-dim", help="PC dimension for reduction", type='integer',default=50)
pca_group$add_argument("--xlow-cutoff", help="x_low_cutoff for search high variable genes, default: 0.125",type='double',default=0.125)
pca_group$add_argument("--xhigh-cutoff", help="x_high_cutoff for search high variable genes, default: 5",type='double',default=5)
pca_group$add_argument("--pc-number", help="PC number for cluster/tSNE/UMAP analysis, default: NULL", type='integer')
cluster_group = parser$add_argument_group("Cluster group", "arguments for cell cluster analysis")
cluster_group$add_argument("--resolution", help="resolution for cluster, greater value and more clusters, default: 0.5", type='double', default=0.5)
cluster_group$add_argument("--annoy-metric", help="Distance metric for annoy. default: euclidean", choices=c("euclidean", "cosine", "manhattan", "hamming"), default="euclidean")
diff_group = parser$add_argument_group("Diff exp group", "arguments for diff analysis or finding markers")
diff_group$add_argument("--min-pct", help="Test genes with min fraction of cells, default: 0.1", type='double',default=0.1)
diff_group$add_argument("--logfc-diff", help="absolute log Foldchange cutoff for differential analysis", type='double',default=0.25)
diff_group$add_argument("--padjust-diff", help="pajust value or power (for roc) cutoff for differential analysis, default: 0.05", type='double',default=0.05)
diff_group$add_argument("--logfc-marker", help="log Foldchange cutoff for marker genes, default: 1", type='double',default=1)
diff_group$add_argument("--padjust-marker", help="pajust value or power (for roc) cutoff for marker genes, default: 0.01", type='double',default=0.01)
diff_group$add_argument("--test-method", help="Denotes which test to use, default: wilcox. If negbinom, poisson, or DESeq2, count data will be used.", choices=c("wilcox","bimod","roc","t","negbinom","poisson","LR","MAST","DESeq2"), default="wilcox")
diff_group$add_argument("--jam-pdf", help="Combine pdf files of cluster-level.", action='store_true')
args = parser$parse_args()
str(args)

input_file = args$input
input_format = args$input_format
prefix = args$prefix
genename_file = args$gene_name
analydir = args$analydir

## qc

# filter paramters
min_cells = args$min_cells
low_feature = args$low_nGene
high_feature = args$high_nGene
high_mt_percent = args$high_pMT
# normalization method
normalize_method = args$normalize_method
genome = args$genome
# mitochondira genes pattern based on reference genome
if(genome == "GRCh38" || genome == "hg19"){
	pattern_MT="^MT"
}
if(genome == "mm10"){
	pattern_MT="^mt-"
}

## pca


selection_method = args$selection_method
pc_dim = args$pc_dim
x_low_cutoff = args$xlow_cutoff
x_high_cutoff = args$xhigh_cutoff
dispersion_cutoff = args$dispersion_cutoff
nfeatures = args$feature_number
pc_number = args$pc_number

## cluster

cluster.resolution = args$resolution
annoy_metric = args$annoy_metric

## diff

logfc.diff = args$logfc_diff
padjust.diff = args$padjust_diff
logfc.marker = args$logfc_marker
padjust.marker = args$padjust_marker
test_method = args$test_method
min_pct = args$min_pc

HVGsdir = file.path(analydir, "Identification_HVGs")
subpopdir = file.path(analydir, "Identification_subpopulation")
Markergenedir = file.path(analydir, "Detection_Markergene")
normalizationdir = file.path(HVGsdir,"Normalization")
pca.dir = file.path(subpopdir, "PCA")
cluster.dir = file.path(subpopdir, "Clustering")
tsne.dir = file.path(subpopdir, "TSNE")
umap.dir = file.path(subpopdir, "UMAP")
dir.create(subpopdir, recursive = TRUE)
dir.create(Markergenedir, recursive = TRUE)
dir.create(normalizationdir, recursive = TRUE)
dir.create(pca.dir, recursive = TRUE)
dir.create(cluster.dir, recursive = TRUE)
dir.create(tsne.dir, recursive = TRUE)
dir.create(umap.dir, recursive = TRUE)

#function
read.file <- function(filename, ...){
        if(grepl(".gz$",filename)){
                return(read.table(gzfile(filename), ...))
        }else{
                return(read.table(filename, ...))
        }
}
Read_genename <- function(genename.file, geneid.col=1, genename.col=2){
        dat = read.file(genename.file, sep="\t", colClasses="character")
        indices = order(dat[,geneid.col])
        GENENAME = data.frame(gid = dat[indices, geneid.col], gname = dat[indices, genename.col], stringsAsFactors=FALSE)
        GENENAME$uname = make.unique(GENENAME$gname)
        return(GENENAME)
}
Read_GeneName <- function(input.file, input.format, genename.file=NULL){
        if(!is.null(genename.file)){
                GENENAME = Read_genename(genename.file)
                return (GENENAME)
        }
        if(input.format == "h5"){
                h5.infile <- hdf5r::H5File$new(filename = input.file, mode = "r")
                GENENAME = data.frame(gid = h5.infile[["matrix/features/id"]][],
                        gname = h5.infile[["matrix/features/name"]][], stringsAsFactors=FALSE)
                GENENAME = GENENAME[order(GENENAME$gid),]
                GENENAME$uname = make.unique(GENENAME$gname)
                return (GENENAME)
        }else if(input.format == "10x" | input.format == "10X"){
                h5.infile <- hdf5r::H5File$new(filename = file.path(input.file,'filtered_feature_bc_matrix.h5'), mode = "r")
                GENENAME = data.frame(gid = h5.infile[["matrix/features/id"]][],
                        gname = h5.infile[["matrix/features/name"]][], stringsAsFactors=FALSE)
                GENENAME = GENENAME[order(GENENAME$gid),]
                GENENAME$uname = make.unique(GENENAME$gname)
                return (GENENAME)
        }else if(input.format == "csvdir"){
                feature.gzfile = file.path(input.file, "features.tsv.gz")
                GENENAME = Read_genename(feature.gzfile)
                return (GENENAME)
        }else{
                return (NULL)
        }
}
## genename colomn names: gene.id, gene.name, unique.name
map.genename <- function(gid, GENENAME=NULL, from=1, to=3){
        if(is.null(GENENAME)){
                return (gid)
        }
        map.vector = GENENAME[,to]
        names(map.vector) = GENENAME[,from]
        x = map.vector[gid]
        indices.na = which(is.na(x))
        x[indices.na] = gid[indices.na]
        names(x) = NULL
        return (x)
}
#
write.xls <- function(df, xls.file, first.colname="Gene", sep="\t", quote=FALSE, ...) {
        dat = data.frame(rownames(df), df, check.names=F)
        colnames(dat)[1] = first.colname
        write.table(dat, file = xls.file, sep=sep, quote=quote, row.names=F, ...)
}
write.xls.genename <- function(df, xls.file, GENENAME=NULL, sep="\t", quote=FALSE, ...) {
        if(is.null(GENENAME)){
                write.xls(df, xls.file, first.colname="Gene", sep="\t", quote=FALSE, ...)
        }else{
                dat = data.frame(Gene = map.genename(rownames(df), GENENAME, from=3, to=1),
                        GeneName = map.genename(rownames(df), GENENAME, from=3, to=2), df, check.names=F)
                write.table(dat, file = xls.file, sep=sep, quote=quote, row.names=F, ...)
        }
}
#
determine_PCnum <- function(pca_fitted, pc.min = 10, pc.max = 20){
        pc.max = min(pc.max, nrow(pca_fitted))
        pc = pc.max
        for (i in 2:pc.max){
                pc.diff = pca_fitted[i-1] - pca_fitted[i]
                if(pc.diff < 0.05){
                        pc = i
                        break
                }
        }
        return(max(pc.min, pc))
}
#
low_res = 70
mid_res = 150
high_res = 300
dual.plot <- function(fig, file.prefix, w=7, res=75){
        pdf(paste(file.prefix,".pdf",sep=""), width = w)
        print(fig)
        dev.off()
        png(paste(file.prefix,".png",sep=""), width = w*res, res = res, type="cairo-png")
        print(fig)
        dev.off()
}
#####################=======================================start to analyze=======================================########################
#####数据再过滤和标准化
#load data
data_10X <- Read10X(data.dir = input_file, gene.column = 2, unique.features = TRUE)

#filter
barcode_matrix <- CreateSeuratObject(counts = data_10X,min.cells = min_cells,min.features = low_feature,project = prefix)
#Calculate mitochondria gene percentage and add this colum as percent.MT into barcode_matrix
barcode_matrix[["percent.MT"]] <- PercentageFeatureSet(barcode_matrix, pattern = pattern_MT)
# filter cells based on gene numbers is lower than high_feature and MT percent is lower than high_mt_percent
barcode_matrix <- subset(x =barcode_matrix, subset=nFeature_RNA < high_feature & percent.MT < high_mt_percent)

cat("normalize ... \n")
feature.names = c("nFeature_RNA","nCount_RNA","percent.MT")
if(normalize_method != "SCT"){
	scale_factor = median(barcode_matrix@meta.data[,feature.names[1]])
	print(paste("scale_factor = ", scale_factor))
        #Normalization method and scale_factor
	barcode_matrix = NormalizeData(object = barcode_matrix, normalization.method = normalize_method, scale.factor = scale_factor)
	barcode_matrix = FindVariableFeatures(object = barcode_matrix,
		selection.method = selection_method,
		mean.function = ExpMean,
		dispersion.function = LogVMR,
		mean.cutoff = c(x_low_cutoff, x_high_cutoff),
		dispersion.cutoff = c(dispersion_cutoff, Inf),
		nfeatures = nfeatures
	)
	sel.features = VariableFeatures(object = barcode_matrix)
        #整理数据ScaleData，使用全部细胞作为后续分析既慢且无必要，features中给出上一步筛选出的表达量变化基因进行后续分析(vars.to.regress)，同时可以设置忽略相关的feature
	barcode_matrix = ScaleData(object = barcode_matrix, features = sel.features, vars.to.regress = feature.names)
        #整理好的数据为后续数据降维和细胞巨雷做准备，去除nUMI和线粒体基因比例，后续降维将不考虑这两个因素，考虑基因标准化后的表达量log（geneUMIcount*10000）
}else{
	barcode_matrix = SCTransform(barcode_matrix, verbose = FALSE, variable.features.n = nfeatures)
	sel.features = VariableFeatures(object = barcode_matrix)
}

# pca 主成分分析，用于多种非随机因素变量进行关联在一起，使用较少的变量来代替之前的变量，线性降维，不能解释特征见的复杂多项式关系，相似数据点不够靠近。 
cat("pca ...\n")
barcode_matrix = RunPCA(object = barcode_matrix, do.print = F, npcs = pc_dim)
## Determine the ‘dimensionality’ 
pc.stdev = data.frame(PC = 1:length(barcode_matrix@reductions$pca@stdev), stdev = barcode_matrix@reductions$pca@stdev)
pc.fit = nls(stdev ~ a*PC^b, data = pc.stdev, start = list(a=10, b= -0.5),trace = T)
if(!is.null(pc_number)){
	pc.num = min(20, pc_number) ## less than 20
	pc.num = max(5, pc.num)    ## greater than 5
}else{
	pc.num = determine_PCnum(fitted(pc.fit))
}
cat(paste("pc.number: ",pc.num,"\n"))
#每个pc随机抽取1%gene，重复100次，得到的p值与均匀分布进行比较
barcode_matrix = JackStraw(barcode_matrix, num.replicate = 100)
barcode_matrix = ScoreJackStraw(barcode_matrix, dims = 1:pc.num)

# cluster  K-nearest neighbor(KNN) graph 进行细胞聚类
# resolution越大得到的细胞分类越多
cat("cluster ...\n")
barcode_matrix = FindNeighbors(barcode_matrix, reduction = "pca", dims = 1:pc.num, annoy.metric = annoy_metric)
barcode_matrix = FindClusters(barcode_matrix, resolution = cluster.resolution)
cluster.ids = levels(barcode_matrix@meta.data$seurat_clusters)

# tSNE / UMAP 非线性降维
# 如果分析目的是鉴定稀有细胞群体，可尝试增加pc的使用个数反复进行下游分析
# 此处只是一种高纬的二维展示，不影响实际聚类结果，umap更能体现全局关系
# T-SNE非线性降维：t- 分布领域嵌入算法 基于pca降维结果，选取pca显著pc进行再次降维，横纵坐标分别代表再次降维后的住成分，不同颜色代表之前聚类中得到的不同分类

cat("tsne / umap ...\n")
barcode_matrix = RunTSNE(object = barcode_matrix, reduction = "pca", dims.use = 1:pc.num, do.fast = TRUE)
barcode_matrix = RunUMAP(object = barcode_matrix, reduction = "pca", dims = 1:pc.num, umap.method = "uwot",metric = "correlation")


##########################################################
###       figures and tables
##########################################################

plot1 = VariableFeaturePlot(barcode_matrix)
plot2 = LabelPoints(plot = plot1, points = sel.features[1:10], repel = TRUE)
#f = CombinePlots(plots = list(plot1, plot2))
dual.plot(plot2, file.path(HVGsdir,paste(prefix,".VariableFeatures.volcano",sep="")), w=12)

cat("read file \n")
GENENAME = Read_GeneName(input.file = input_file, input.format = input_format, genename.file = genename_file)
has_genename = TRUE
if(is.null(GENENAME)){
        has_genename = FALSE
}

#counts.dat = round(as_label(barcode_matrix[["RNA"]]@counts,"matrix"),5)
write.xls.genename(round(barcode_matrix[["RNA"]]@counts,5),
	file.path(normalizationdir, paste(prefix,".UMI_count.matrix.xls",sep="")),
	GENENAME)
write.xls.genename(round(barcode_matrix[["RNA"]]@data,5), 
	file.path(normalizationdir, paste(prefix,".normalize_expression.matrix.xls",sep="")), 
	GENENAME)

HVGgene = data.frame(Gene=sel.features, Genename=map.genename(sel.features,GENENAME,3,1))
write.table(HVGgene,file=file.path(HVGsdir,paste(prefix,".VariableFeatures.genes.txt",sep="")),quote=F,sep="\t",row.names=F)

## pca out files

#查看特定pc中的基因和基因对应的贡献度
f = VizDimLoadings(barcode_matrix, dims = 1:2, reduction = "pca")
dual.plot(f, file.path(pca.dir,paste(prefix,".PCA.topGene",sep="")), w=8, res=mid_res)

# DimPlot()绘制pca散点图，一般只绘制pc1和pc2

#根据pca结果绘制基因表达量热图，dims可以指定绘制那个pc，cells设置使用细胞的个数，可以直观查看主成分的代表性，
#可用于辅助判断后续聚类和降维使用的pc个数，一般选取1～15个，细胞类群较多可以尝试使用更多pc
f = DimHeatmap(barcode_matrix, dims = 1:15, cells = 500, balanced = TRUE, fast=FALSE, raster=FALSE)
dual.plot(f, file.path(pca.dir,paste(prefix,".PCA.heatmap",sep="")), w=12, res=mid_res)

p = ggplot(pc.stdev, aes(PC,stdev)) + geom_point(size=3) + geom_line(aes(PC,fitted(pc.fit)),col='red')
p = p+xlab("PC")+ylab("Standard Deviation of PC")
dual.plot(p, file.path(pca.dir,paste(prefix,".PCA.sdev_fitted",sep="")), w=8, res=mid_res)


# 通过判断图中的拐点来判断选择前多少个pc继续进行下游分析，如果分析目的是鉴定稀有细胞群体，可尝试增加pc的使用个数反复进行下游分析
f = ElbowPlot(barcode_matrix)
dual.plot(f, file.path(pca.dir,paste(prefix,".PCA.ElbowPlot",sep="")), w=8, res=mid_res)


# 每个pc随机抽取1%基因，重复100次，得到的p值与均匀分布进行比较。比较显著的pc会聚集在p值比较小的部分
f = JackStrawPlot(barcode_matrix, dim = 1:15)
dual.plot(f, file.path(pca.dir,paste(prefix,".PCA.JackStrawPlot",sep="")), w=8, res=mid_res)

#########   non-linear dimensional reduction, tSNE / UMAP

write.xls(barcode_matrix@reductions$tsne@cell.embeddings, file.path(tsne.dir,paste(prefix,".cluster.tSNE.xls",sep="")), first.colname="barcode")
tsneplot = DimPlot(barcode_matrix, reduction = "tsne",label = TRUE)
dual.plot(tsneplot, file.path(tsne.dir,paste(prefix,".cluster.tSNE",sep="")), w=8, res=mid_res)

write.xls(barcode_matrix@reductions$umap@cell.embeddings, file.path(umap.dir,paste(prefix,".cluster.UMAP.xls",sep="")), first.colname="barcode")
umapplot = DimPlot(barcode_matrix, reduction = "umap",label = TRUE)
dual.plot(umapplot, file.path(umap.dir,paste(prefix,".cluster.UMAP",sep="")), w=8, res=mid_res)

#### cluster
cluster.res = barcode_matrix@meta.data[,c("orig.ident","seurat_clusters")]
write.xls(cluster.res, file.path(cluster.dir, paste(prefix,".clusters_info.xls",sep="")), first.colname="barcode")

########  biomarkers / diff analysis ####################
# FindAllMarkers可以自动得到每个聚类和其他聚类相比的marker gene
# FindMarkers也可以用ident.1, ident.2指定比较的聚类，min.pct设置基因至少要出现在两个做比较聚类的**%的细胞中
# 差异检验的方式：test.use参数指定

diff.exp = FindAllMarkers(
	object = barcode_matrix, 
	only.pos = F, 
	min.pct = min_pct,
	test.use = test_method,
	logfc.threshold = 0.1,
	return.thresh = 1
)
if(has_genename){
	diff.exp$gene_id = map.genename(diff.exp$gene, GENENAME, 3, 1)
	diff.exp$gene = map.genename(diff.exp$gene_id, GENENAME, 1, 2)
}

diff.exp = diff.exp[,c(8,7,6,2,1,5,3,4)]
diff.sig.exp = subset(diff.exp, abs(avg_logFC)>logfc.diff & p_val_adj<padjust.diff)
marker.exp = subset(diff.exp, avg_logFC>logfc.marker & p_val_adj<padjust.marker)
# top6 平均差异倍数最大的6个marker gene
top6.markers <- diff.exp %>% group_by(cluster) %>% top_n(n = 6, wt = avg_logFC)
# top10 平均差异倍数最大的10个marker gene
top10.markers <- diff.exp %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

diff.file = file.path(Markergenedir,paste(prefix,".cluster.diffexp.xls",sep=""))
diff.sig.file = file.path(Markergenedir,paste(prefix,".cluster.diffexp.significant.xls",sep=""))
marker.file = file.path(Markergenedir,paste(prefix,".cluster.markers.xls",sep=""))
marker.top10.file = file.path(Markergenedir,paste(prefix,".cluster.top10.markers.xls",sep=""))
write.table(diff.exp, file=diff.file, quote=F, sep="\t", row.names=F)
write.table(diff.sig.exp, file=diff.sig.file, quote=F, sep="\t", row.names=F)
write.table(marker.exp, file=marker.file, quote=F, sep="\t", row.names=F)
#write.table(top10.markers, file=marker.top10.file, quote=F, sep="\t", row.names=F)

# marker基因在各个细胞中的表达热图，此处每个cluster只展示平均差异倍数最大的10个marker基因
f = DoHeatmap(barcode_matrix, features = top10.markers$gene)+scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill") #+ NoLegend()
dual.plot(f, file.path(Markergenedir,paste(prefix,".cluster.top10.genename.heatmap",sep="")), w=12, res=high_res)
f = DoHeatmap(barcode_matrix, features = top10.markers$gene, raster=FALSE)+scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill") #+ NoLegend()
dual.plot(f, file.path(Markergenedir,paste(prefix,".cluster.top10.heatmap",sep="")), w=12, res=high_res)

AverageExp = AverageExpression(barcode_matrix, features=unique(top10.markers$Gene_name))
corelation = corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
f = pheatmap(corelation$r)
dual.plot(f, file.path(Markergenedir,paste(prefix,".cluster.correlation",sep="")), res=high_res)

jam.pdfs = NULL
for(i in cluster.ids){
	## violin plot
	genes.plot = top6.markers$gene[top6.markers$cluster==i]
	v = VlnPlot(barcode_matrix, genes.plot, ncol =2, pt.size = 0) + xlab("Cluster") + ylab("log(UMI)")
	dual.plot(v, file.path(Markergenedir,paste(prefix,".cluster",i,".marker.violin",sep="")), w=10, res=mid_res)
	## feature plot
	f1 = FeaturePlot(object = barcode_matrix, genes.plot, cols = c("grey", "blue"), reduction = "tsne")
	dual.plot(f1, file.path(Markergenedir,paste(prefix,".cluster",i,".marker.tsne",sep="")), res=mid_res)
	f2 = FeaturePlot(object = barcode_matrix, genes.plot, cols = c("grey", "blue"), reduction = "umap")
	dual.plot(f2, file.path(Markergenedir,paste(prefix,".cluster",i,".marker.umap",sep="")), res=mid_res)
	jam.pdfs = c(jam.pdfs, file.path(Markergenedir,paste(prefix,".cluster",i,".marker.tsne.pdf",sep="")),
		file.path(Markergenedir,paste(prefix,".cluster",i,".marker.umap.pdf",sep="")),
		file.path(Markergenedir,paste(prefix,".cluster",i,".marker.violin.pdf",sep="")))
}
if(args$jam_pdf){
	jamPDF(jam.pdfs,
		out.file = file.path(Markergenedir,paste(prefix,".cluster.markers",sep="")),
		layout = '2x2',
		delete.original = FALSE,
		ignore.stderr = FALSE)
}
