##########Cell type annotation using SingleR
#Author: Dongfang Hu
#Date: 20200621
#####################################


library(argparse)
library(Seurat)
library(SingleR)
library(cowplot)
library(ggplot2)

parser = ArgumentParser(description="SingleR analysis for sc-RNAseq")
parser$add_argument("--input", help="input files or directories. Required.", nargs='+', required=T)
parser$add_argument("--prefix", help="prefix for output files. Required.", required=T)
parser$add_argument("--analydir", help="analysis directory. Required.", required=T)
parser$add_argument("--species", help="analysis species. Required.", required=T)
args = parser$parse_args()
str(args)

input = args$input
prefix = args$prefix
analydir = args$analydir
species = args$species
outdir = file.path(analydir,prefix)
dir.create(outdir, recursive = TRUE)

if(species == "hsa"){
	refspecies = "Human"
	reflist = list(hpca=hpca, blueprint_encode=blueprint_encode)
	reflistname = c("hpca","blueprint_encode")
} else if (species == "mmu"){
	refspecies = "Mouse"
	reflist = list(immgen=immgen,mouse.rnaseq=mouse.rnaseq)
	reflistname = c("immgen","mouse.rnaseq")
} else {
	stop("your species is not hsa or mmu, please prepare refdata")
}
#function
dual.plot <- function(fig, file.prefix){
        pdf(paste(file.prefix,".pdf",sep=""))
        print(fig)
        dev.off()
        png(paste(file.prefix,".png",sep=""), type="cairo-png")
        print(fig)
        dev.off()
}
dual.plot1 <- function(fig, file.prefix){
        pdf(paste(file.prefix,".pdf",sep=""))
        print(fig)
        dev.off()
}

findCells <- function(obj, column, values, name=NULL) {
  ## Given a Seurat OBJ, return a list with the names of the cells where
  ## the specified meta.data COLUMN equals any of the strings specified
  ## in VALUES (both must be characters or factors). Name of the list member
  ## must be specified using the NAME argument if length(values)>1
  stopifnot(is(obj, "Seurat"))
  stopifnot(is.character(column))
  stopifnot(column %in% names(obj@meta.data))
  col <- obj@meta.data[[column]]
  stopifnot(is.character(col) || is.factor(col))
  values <- unique(values)
  stopifnot(is.character(values) || is.factor(values))
  if (length(values)>1 && is.null(name))
    stop("findCells: specify a name to be used for the selection")
  if(is.null(name))
    name <- values
  stopifnot(is.character(name))
  rem <- setdiff(c(values), col)
  if(length(rem)>0)stop("findCells: requested value(s) never occurs in this column: ", rem)
  l <- list(colnames(obj)[ col %in% values ])
  names(l) <- name
  l
}

#load Seuratobject
Seuratobject <- readRDS(file=input)
counts <- GetAssayData(Seuratobject)

print("Star singler")
singler <- CreateSinglerObject(counts=counts,
  project.name="10X scRNAseq", # choose
  min.genes = 0, # ignore cells with fewer than 0 genes
  technology = "10X",
  species = refspecies, #choose Human or Mouse
  citation = "",
  ref.list = reflist, # hsa: hpca and blueprint_encode; mmu: immgen and mouse.rnaseq
  normalize.gene.length = FALSE,  # needed for full-length platforms (e.g. smartseq)
  variable.genes = "de",  # 'de' uses pairwise difference between the cell types
  fine.tune = FALSE, # TRUE would take very long
  reduce.file.size = TRUE, # leave out less-often used fields 
  do.signatures = FALSE,
  do.main.types = TRUE,
  numCores = SingleR.numCores)

write.table(singler$singler[[1]]$SingleR.single.main, file = file.path(outdir,paste(prefix,"_",reflistname[1],"_maintype.xls",sep="")), quote=F,sep="\t",row.names=F)
write.table(singler$singler[[2]]$SingleR.single.main, file = file.path(outdir,paste(prefix,".",reflistname[2],"_maintype.xls",sep="")), quote=F, sep="\t",row.names=F)

#only show in screen
for (ref.set in names(singler$singler) ) {
  types <- singler$singler[[ref.set]]$SingleR.single.main$labels[,1]
  cat("==== ", ref.set, ": ====\n")
  show(sort(table(types), decreasing=TRUE))
}

#add singler cell lable to Seuratobject
for (ref.set in names(singler$singler) ) {
  types <- singler$singler[[ref.set]]$SingleR.single.main$labels[,1]
  Seuratobject <- AddMetaData(Seuratobject,
                         metadata=types,
                         col.name=paste0(ref.set,"_type" ) )
}

write.table(Seuratobject@meta.data,file=file.path(outdir, paste(prefix,".Celltype_annotation.xls",sep='')),sep='\t',quote=FALSE,row.names=T)

print("plot TSEN/UMAP...")
#plot
p1 <- DimPlot(Seuratobject, group.by=paste(reflistname[1],'_type',sep=''), reduction="tsne")
p2 <- DimPlot(Seuratobject, group.by=paste(reflistname[2],'_type',sep=''), reduction="tsne")
p <- plot_grid(p1, p2, nrow=2, ncol=1, labels=reflistname)
dual.plot(p, file.path(outdir,paste(prefix,".Celltype_TSNE",sep="")))
p3 <- DimPlot(Seuratobject, group.by=paste(reflistname[1],'_type',sep=""), reduction="umap")
p4 <- DimPlot(Seuratobject, group.by=paste(reflistname[2],'_type',sep=""), reduction="umap")
f <- plot_grid(p3, p4, nrow=2, ncol=1, labels=reflistname)
dual.plot(f, file.path(outdir,paste(prefix,".Celltype_UMAP",sep="")))

#B cell plot
print("plot B cell...")
p5 <- DimPlot(Seuratobject, group.by=paste(reflistname[1],'_type',sep=""), reduction="tsne", cells.highlight=findCells(Seuratobject, paste(reflistname[1],'_type',sep=""),
                 c('B_cell', 'Pre-B_cell_CD34-', 'Pro-B_cell_CD34+'),
                 name="B-like"))
p6 <- DimPlot(Seuratobject, group.by=paste(reflistname[2],'_type',sep=""), reduction="tsne", cells.highlight=findCells(Seuratobject, paste(reflistname[2],'_type',sep=""),
                 c('B-cells'),
                 name="B-cells"))
f <- plot_grid(p5, p6, nrow=2, ncol=1, labels=paste(reflistname, "all B cell"))
dual.plot(f, file.path(outdir,paste(prefix,".TSNE_Bcells",sep="")))

p7 <- DimPlot(Seuratobject, group.by=paste(reflistname[1],'_type',sep=""), reduction="umap", cells.highlight=findCells(Seuratobject, paste(reflistname[1],'_type',sep=""),
                 c('B_cell', 'Pre-B_cell_CD34-', 'Pro-B_cell_CD34+'),
                 name="B-like"))
p8 <- DimPlot(Seuratobject, group.by=paste(reflistname[2],'_type',sep=""), reduction="umap", cells.highlight=findCells(Seuratobject, paste(reflistname[2],'_type',sep=""),
                 c('B-cells'),
                 name="B-cells"))
f <- plot_grid(p7, p8, nrow=2, ncol=1, labels=paste(reflistname, "all B cell"))
dual.plot(f, file.path(outdir,paste(prefix,".UMAP_Bcells",sep="")))

#T cell plot
p9 <- DimPlot(Seuratobject, group.by=paste(reflistname[1],'_type',sep=""),reduction="tsne", cells.highlight=findCells(Seuratobject, paste(reflistname[1],'_type',sep=""),
		c('T_cells'),
		name="T_cells"))
p10 <- DimPlot(Seuratobject, group.by=paste(reflistname[2],'_type',sep=""), reduction="tsne", cells.highlight=findCells(Seuratobject, paste(reflistname[2],'_type',sep=""),
		c('CD4+ T-cells','CD8+ T-cells'),
                name="T_cells"))
f <- plot_grid(p9, p10, nrow=2, ncol=1, labels=paste(reflistname, "all T cell"))
dual.plot(f, file.path(outdir,paste(prefix,".TSNE_Tcells",sep="")))

p11 <- DimPlot(Seuratobject, group.by=paste(reflistname[1],'_type',sep=""), reduction="umap", cells.highlight=findCells(Seuratobject, paste(reflistname[1],'_type',sep=""),
                c('T_cells'),
                name="T_cells"))
p12 <- DimPlot(Seuratobject, group.by=paste(reflistname[2],'_type',sep=""), reduction="umap", cells.highlight=findCells(Seuratobject, paste(reflistname[2],'_type',sep=""),
                c('CD4+ T-cells','CD8+ T-cells'),
                name="T_cells"))
f <- plot_grid(p11, p12, nrow=2, ncol=1, labels=paste(reflistname, "all T cell"))
dual.plot(f, file.path(outdir,paste(prefix,".UMAP_Tcells",sep="")))

#####
singler$seurat = Seuratobject
singler$meta.data$orig.ident = Seuratobject@meta.data$orig.ident
singler$meta.data$xy = Seuratobject@reductions$tsne@cell.embeddings
singler$meta.data$xy = Seuratobject@reductions$umap@cell.embeddings
singler$meta.data$clusters = Seuratobject@active.ident
save(singler,file=file.path(outdir,paste(prefix,'.SingleRobject.RData',sep="")))
write.table(singler$meta.data,file=file.path(outdir,paste(prefix,".SingleRobject_metadata.xls",sep="")),quote=F,sep="\t",row.names=T)

Table1 = table(singler$meta.data$orig.ident,singler$singler[[1]]$SingleR.single.main$labels)
write.table(Table1, file=file.path(outdir,paste(prefix,".",reflistname[1],"_number.xls",sep="")), quote=F,sep="\t",row.names=T)
Table2 <- table(singler$meta.data$orig.ident,singler$singler[[2]]$SingleR.single.main$labels)
write.table(Table2, file=file.path(outdir,paste(prefix,".",reflistname[2],"_number.xls",sep="")), quote=F,sep="\t",row.names=T)


f <- SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, normalize = T, top.n = 50, clusters = singler$meta.data$clusters,order.by.clusters=T)
dual.plot(f, file.path(outdir,paste(prefix,".",reflistname[1],".main.heatmap_top50",sep="")))
f <- SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, normalize = T,clusters = singler$meta.data$clusters,order.by.clusters=T)
dual.plot1(f, file.path(outdir,paste(prefix,".",reflistname[1],".main.heatmap_detail",sep="")))
f <- SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, normalize = T, top.n = 50, clusters = singler$meta.data$clusters,order.by.clusters=T)
dual.plot(f, file.path(outdir,paste(prefix,".",reflistname[1],".heatmap_top50",sep="")))

f <- SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single.main, normalize = T, top.n = 50, clusters = singler$meta.data$clusters,order.by.clusters=T)
dual.plot(f, file.path(outdir,paste(prefix,".",reflistname[2],".main.heatmap_top50",sep="")))
f <- SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single.main, normalize = T, clusters = singler$meta.data$clusters,order.by.clusters=T)
dual.plot1(f, file.path(outdir,paste(prefix,".",reflistname[2],".main.heatmap_detail",sep="")))
f <- SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single, normalize = T, top.n = 50, clusters = singler$meta.data$clusters,order.by.clusters=T)
dual.plot(f, file.path(outdir,paste(prefix,".",reflistname[2],".heatmap_top50",sep="")))
