###
# Author Dongfang Hu
# This script is designed to pipeline Seurat object output into Monocle2
# It will start from reading in a Seurat object with the count sparse matrix, UMAP coordinates, and cluster information
# Date 20200701
####



#######################单细胞拟时间分析####################
#scRNA-seq数据的表达矩阵之后大多数后续分析的重点是：确定组织或者癌症中细胞亚群的类型或状态，或研究动态变化过程，例如细胞分化，细胞周期或刺激反应
# 使用聚类算法将细胞分为不同的细胞类型或状态
# 通过轨迹推断方法沿伪时间轴对细胞进行排序
#细胞在不同的状态下都有特异行驶功能的基因
#通过关键基因的表达，对细胞进行排序，进而模拟出随着时间的进行，细胞的发育状态的改变
#轨迹推断的前提
#对目标细胞亚群已经有清晰的鉴定和划分
#确保基因集和降维符合预期
#预期其中若干细胞亚群存在谱系上的分化继承关系，并期望通过分析确定它们的分化过程
#轨迹推断的目的：推导若干个可能存在分化/演化继承关系的细胞亚群间最可能的分化路径
#找出驱动细胞亚群分化的关键基因


#############################monocle原理介绍################
#发育轨迹分析，也称为拟时间分析，一般使用monocle2软件进行分析，monocle2软件是一个擅长进行发育轨迹分析的综合工具，它可以进行质控，聚类，差异分析，关键节点基因筛选等功能
#通过拟时间分析，我们可以得到细胞的发育轨迹，以及发育过程中的关键基因的表达变化，找到发育过程中起关键作用的基因。
#如何选取基因集：an unsupervised procedure "dpFeature", selecting the genes that are differentially expressed between clusters of cells identified with t-distributed stochastic neighbor embedding (tSNE) dimension reduction followed by density-peak clustering.
#如何排序: Monocle ICA (v1) Monocle 1 默认采用ICA算法进行数据降维，ICA算法使用MST（Minimum spamming tree）获得细胞间分化轨迹的最小路径 缺点1，假设细胞间是无关系的，2假设细胞呈现高斯分布
        # Monocle DDRTree (v2) Monocle2默认采用DDRTree算法进行数据降维（更适用于细胞轨迹分析）
                #筛选细胞类cluster的所有差异表达基因，降维并构建最小生成树，对单细胞进行在高维和低纬空间搜索最优排序，拟合最佳细胞发育或者分化你是轨迹曲线
        #Monocle UMAP （v3），处理的细胞数增加，支持UMAP推断发育轨迹， 软件安装比较麻烦
#如何确定分支

###########################monocle分析流程###################
#scRNAseq dataset
#Pre-process data: Normalize Remove batch effects
#Non-linear dimensionality reduction
#cluster cells
#compare clusters Identify top markers, Targeted contrasts
#Trajectory analysis





library(Seurat)
library(monocle)
library(argparse)
parser = ArgumentParser(description="Cell trajectory analysis")
parser$add_argument("--input", help="input files or directories. Required.", required=T)
parser$add_argument("--outdir", help="analysis directory. Required.", required=T)
parser$add_argument("--prefix", help="sample or group",required=T)
args = parser$parse_args()
str(args)

inputfile = args$input
outdir = args$outdir
prefix = args$prefix

dir.create(outdir,showWarnings = FALSE,recursive = TRUE)

#function 绘制两种格式图片功能pdf和png
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

#Extract data, phenotype data, and feature data from the SeuratObject 从seurat2monocle构建
data <- as(as.matrix(seurat@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds 创建monocle2 CDS 对象 
print("Construct monocle newCellDataSet...")
cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              #lowerDetectionLimit = 0.5, UMI， 这一步在之前的suerat分析中已经去去除所以不需要再筛选
                              expressionFamily = uninormal())# since I have already normalized, thresholded and scalled in Suerat v3.


#Run ordering algorithm 
var_genes <- seurat[["RNA"]]@var.features
ordering_genes <- var_genes
cds <- setOrderingFilter(cds, ordering_genes)
print(dim(exprs(cds)))

## reduce dimension - do not normalize or include pseudo count. Use monocle scaling 降维
cds <- reduceDimension(cds,norm_method="none", 
                        reduction_method="DDRTree",
                        max_components=4,
                        scaling=TRUE,
                        verbose=TRUE,
                        pseudo_expr=0)

# First decide what you want to color your cells by
print(head(pData(cds)))

## order cells change colors and theta to match your plot 每个细胞按照轨迹进行排序
cds <- orderCells(cds)
write.table(pData(cds), file = file.path(outdir,paste(prefix,".Trajectory_BarcodeMatrix.xls",sep="")), sep='\t', quote=F, row.names=T)

#plot 绘制轨迹 color_by = "State", color_by = "seurat_clusters" + facet_wrap(~XXX, nrow = 1), color_by = "Pseudotime"
p <- plot_cell_trajectory(cds, 
                     color_by = "seurat_clusters",
                     theta = -15,
                     show_tree = TRUE)
dual.plot2(p, file.path(outdir,paste(prefix,".TrajectoryPlot_by_clusters",sep="")),w=8, h=5,res=150)

p <- plot_cell_trajectory(cds, color_by = "Pseudotime",theta = -15,show_tree = TRUE)
dual.plot2(p, file.path(outdir,paste(prefix,".TrajectoryPlot_by_Pseudotime",sep="")),w=8, h=5,res=150)



#查看核心基因轨迹, 也可以查看核心基因表达热图plot_pseudotime_heatmap
my_genes <- head(row.names(subset(fData(cds))))
cds_subset <- cds[my_genes, ]
f <- plot_genes_in_pseudotime(cds_subset, color_by = "orig.ident")
dual.plot2(f, file.path(outdir,paste(prefix,"Plot_genes_in_Pseudotime",sep="")),w=8, h=5,res=150)
#也可以查看核心基因表达热图plot_pseudotime_heatmap

#BEAM分析 不同分支的基因热图 BEAM plot_genes_branched_heatmap plot_genes_branched_pseudotime