library(Seurat)
library(ggpubr)
library(dplyr)
library(clustree)
library(monocle)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(corrplot)
library(org.Dm.eg.db)
library(topGO)
library(biomaRt)
library(scater)
library(velocyto.R)
library(SeuratWrappers)
library(DropletUtils)
library(WGCNA)
theme_set(theme_cowplot())


ldat <- ReadVelocity(file = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/Ctrl_FC/Ctrl_Ovary.loom")
pbmc <- as.Seurat(x = ldat)

ggplot(pbmc@meta.data, aes(nFeature_spliced)) + geom_histogram(bins=100, color = "black", fill = "white") + xlab("Sum of expression") + ggtitle("Total expression histogram before normalization")
ggplot(pbmc@meta.data, aes(nFeature_spliced)) + geom_density(alpha = 0.2, color = "red", fill = "red") + xlab("Sum of expression") + ggtitle("Total expression density before normalization")
ggplot(pbmc@meta.data, aes(nCount_spliced)) + geom_histogram(bins=100, color = "black", fill = "white") + xlab("Total nCount_spliced") + ggtitle("Total UMIs histogram before normalization")
ggplot(pbmc@meta.data, aes(nCount_spliced)) + geom_density(alpha = 0.2, color = "red", fill = "red") + xlab("Total nCount_spliced") + ggtitle("Total UMIs density before normalization")
ggplot(pbmc@meta.data, aes(nFeature_spliced, nCount_spliced)) + geom_point() + geom_smooth(method=lm) + stat_cor(method = "pearson", label.x = 4000, label.y = 3000) + xlab("nFeature_spliced") + ylab("nCount_spliced")

table(pbmc@meta.data$orig.ident)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "-m")
pbmc[["percent.rp"]] <- PercentageFeatureSet(pbmc, pattern = "^Rp")
VlnPlot(pbmc, features = c("nFeature_spliced", "nCount_spliced", "percent.mt", "percent.rp"), ncol = 2)
p1 <- FeatureScatter(pbmc, feature1 = "nCount_unspliced", feature2 = "nCount_spliced")
p2 <- FeatureScatter(pbmc, feature1 = "nFeature_unspliced", feature2 = "nFeature_spliced")
p3 <- FeatureScatter(pbmc, feature1 = "nCount_spliced", feature2 = "nFeature_spliced")
p4 <- FeatureScatter(pbmc, feature1 = "nCount_unspliced", feature2 = "nFeature_unspliced")
CombinePlots(plots = list(p1, p2, p3, p4))

pbmc <- subset(pbmc, subset = nFeature_spliced > 775 & nFeature_spliced < 2200 & nCount_spliced < 18000 & nCount_spliced > 2000 & percent.mt < 5)
pbmc <- NormalizeData(pbmc)

#Assign Limits
ggplot(pbmc@meta.data, aes(nFeature_spliced)) + geom_density(alpha = 0.2, color = "red", fill = "red") + xlim(0,4000) + xlab("Sum of expression") + ggtitle("Total expression density before normalization")
ggplot(pbmc@meta.data, aes(nCount_spliced)) + geom_density(alpha = 0.2, color = "red", fill = "red") + xlim(0,45000) + xlab("Total nCount_spliced") + ggtitle("Total UMIs density before normalization")
ggplot(pbmc@meta.data, aes(nFeature_spliced, nCount_spliced)) + geom_point() + geom_smooth(method=lm) + stat_cor(method = "pearson", label.x = 4000, label.y = 3000) + xlab("nFeature_spliced") + ylab("nCount_spliced") + geom_vline(xintercept = c(500,1000,2000,3000,4000), linetype="dashed", color="grey") + geom_hline(yintercept = c(1000,45000), linetype="dashed", color="grey")
table(pbmc@meta.data$orig.ident)
ggplot(pbmc@meta.data, aes(pbmc@meta.data$orig.ident,nFeature_spliced)) + geom_violin(trim = FALSE, adjust = 1, aes(fill=pbmc@meta.data$orig.ident)) + geom_hline(yintercept = c(1000,2000,3000)) +NoLegend()
#ASSIGN_LIMS_TO_PRIM_ImR
ggplot(pbmc@meta.data, aes(nFeature_spliced)) + geom_density(alpha = 0.2, color = "red", fill = "red") + xlim(0,4000) + xlab("Sum of expression") + ggtitle("Total expression density before normalization") + geom_vline(xintercept = c(200,1000,1400,2000,2350,2750,3150,3600,4000), linetype="dashed", color="grey")
ggplot(pbmc@meta.data, aes(nFeature_spliced, nCount_spliced)) + geom_point() + geom_smooth(method=lm) + stat_cor(method = "pearson", label.x = 2500, label.y = 3000) + xlab("nFeature_spliced") + ylab("nCount_spliced") + geom_vline(xintercept = c(200,1000,1400,2000,2350,2750,3150,3600,4000), linetype="dashed", color="grey") + geom_hline(yintercept = c(1000,45000), linetype="dashed", color="grey")
table(pbmc@meta.data$orig.ident)
ggplot(pbmc@meta.data, aes(pbmc@meta.data$orig.ident,nFeature_spliced)) + geom_violin(trim = FALSE, adjust = 0.5, aes(fill=pbmc@meta.data$orig.ident)) + geom_hline(yintercept = c(1000,1400,2000,2350,2750,3150,3600))

ggplot(pbmc@meta.data, aes(nFeature_spliced)) + geom_density(alpha = 0.2, color = "red", fill = "red") + xlab("Sum of expression") + ggtitle("Total expression density after normalization")
ggplot(pbmc@meta.data, aes(nCount_spliced)) + geom_density(alpha = 0.2, color = "red", fill = "red") + xlab("Total nFeature_spliced") + ggtitle("Total UMIs density after normalization") #+ xlim(650,2500)
table(pbmc@meta.data$orig.ident)
VlnPlot(pbmc, features = c("nFeature_spliced", "nFeature_spliced", "percent.mt"), pt.size = 0.8, ncol = 3)
ggplot(pbmc@meta.data, aes(nFeature_spliced, nCount_spliced)) + geom_point() + geom_smooth(method=lm) + stat_cor(method = "pearson", label.x = 1500, label.y = 3000) + xlab("nFeature_spliced") + ylab("nCount_spliced")
table(pbmc@meta.data$orig.ident)

cc.genes <- readLines(con = "/home/dchatterjee/cell_cycle_genes.txt")
g2m.genes <- cc.genes[1:68]
s.genes <- cc.genes[69:124]
pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pbmc$CC.Difference <- pbmc$S.Score - pbmc$G2M.Score

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#pbmc <- FindVariableFeatures(pbmc, selection.method = "disp", nfeatures = 2000) #THIS IS THE v2.3.4 METHOD
#pbmc <- FindVariableFeatures(pbmc, selection.method = "mvp", nfeatures = 2000) #ANOTHER OPTION
top30 <- head(VariableFeatures(pbmc), 30)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

pbmc <- SCTransform(pbmc, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_spliced"), assay = "spliced", verbose = TRUE)

pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
DimPlot(pbmc, reduction = "pca", group.by = "orig.ident")
VizDimLoadings(pbmc, dims = 1:9, reduction = "pca")
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#pbmc <- JackStraw(pbmc, num.replicate = 50, dims = 100)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:100)
#JackStrawPlot(pbmc, dims = 1:100)
ElbowPlot(pbmc, ndims = 100)

pbmc <- FindNeighbors(pbmc, dims = 1:75)
pbmc <- FindClusters(pbmc, resolution = 0.5)
table(pbmc@active.ident,pbmc@meta.data$orig.ident)
pbmc <- RunUMAP(pbmc, dims = 1:30, n.neighbors = 20, min.dist = 0.35, spread = 0.5)
DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#Merge Clusters
Idents(object = pbmc, cells = WhichCells(pbmc, idents = "22")) <- "21"

cluster.markers <- FindMarkers(pbmc, ident.1 = 5, logfc.threshold = 0.25, test.use = "DESeq2", only.pos = TRUE)
head(cluster.markers, n = 15)

FeaturePlot(pbmc, features = c("Ilp6", "Past1", "ct", "peb"), pt.size = 0.5, cols = c("lightgrey", "red"))

pbmc <- RunUMAP(pbmc, dims = 1:75, n.neighbors = 150, n.epochs = 1000, min.dist = 0.45, repulsion.strength = 0.8, spread = 1.2)
DimPlot(pbmc, reduction = "umap", pt.size = 1, label = TRUE) #2D
pbmc <- RunUMAP(pbmc, dims = 1:75, n.neighbors = 150, n.components = 3, n.epochs = 1000, min.dist = 0.45, repulsion.strength = 0.8, spread = 1.2)
DimPlot(pbmc, reduction = "umap", pt.size = 1, label = TRUE) #3D

bm <- RunVelocity(object = pbmc, deltaT = 1, kCells = 20, min.nmat.emat.correlation = 0.2, min.nmat.emat.slope = 0.2, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)

x <- show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)

emat <- pbmc[["spliced"]]
nmat <- pbmc[["unspliced"]]

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = pbmc)))
names(x = ident.colors) <- levels(x = pbmc)
cell.colors <- ident.colors[Idents(object = pbmc)]
names(x = cell.colors) <- colnames(x = pbmc)

gene <- "ct"
par(mfrow=c(2,2))
gene.relative.velocity.estimates(nmat, emat, deltaT=1,kCells = 100,kGenes=20,fit.quantile=0.2,cell.emb=Embeddings(object = pbmc, reduction = "umap"),cell.colors=cell.colors,show.gene=gene,do.par=T)

##LglIR and CTRL

pbmc <- FindNeighbors(pbmc, dims = 1:50)
pbmc <- FindClusters(pbmc, resolution = 1.5)
table(pbmc@active.ident,pbmc@meta.data$orig.ident)
pbmc <- RunUMAP(pbmc, dims = 1:50, n.neighbors = 25, min.dist = 0.3, spread = 0.65)
DimPlot(pbmc, reduction = "umap", pt.size = 1, label = TRUE)

#TO GET NORMALIZED COUNT FOR A GENE FROM ALL CELLS IN PBMC
library(lava)
DefaultAssay(pbmc) <- "RNA" #CHOOSE ASSAY
dataset <- t(as.matrix(GetAssayData(object = pbmc, slot = "counts"))) #RAW COUNT
dataset <- t(as.matrix(GetAssayData(object = pbmc, slot = "data"))) #NORMALIZED COUNT
max(as.matrix(Grep(as.data.frame(dataset), "Ets21C", subset = TRUE, ignore.case = TRUE))) #CHOOSE MIN/MEDIAN/MAX

#TO GET NORMALIZED COUNTS FOR A GENE FROM CELLS IN A PARTICULAR CLUSTER IN PBMC
subset <- subset(pbmc, idents = c("23"))
DefaultAssay(subset) <- "RNA"
dataset <- t(as.matrix(GetAssayData(object = subset, slot = "data")))
mean(as.matrix(Grep(as.data.frame(dataset), "Ets21C", subset = TRUE, ignore.case = TRUE)))

##WGCNA

options(stringsAsFactors = F)
datExpr <- t(as.matrix(GetAssayData(object = pbmc, slot = "data")))[,VariableFeatures(pbmc)]


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net <- blockwiseModules(datExpr, power = 26, deepSplit = 3,
						corType = "bicor",
						networkType = "unsigned", minModuleSize = 20,
						reassignThreshold = 0, mergeCutHeight = 0.15,
						numericLabels = TRUE, pamRespectsDendro = FALSE,
						saveTOMs = TRUE,
						saveTOMFileBase = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/LglKD+CtrlWT_FC/WGCNA/ALRA_beta26_deepSplit3_unsigned_minMod20_cutHeight015_TOM",
						verbose = 3)

table(net$colors)

mergedColors = net$colors

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

f <- function(module){
    eigengene <- unlist(net$MEs[paste0("ME", module)])
    means <- tapply(eigengene, Idents(pbmc), mean, na.rm = T)
    return(means)
}

colnames(datExpr)[net$colors == "blue"]

modules <- c("blue", "brown", "green", "grey", "turquoise", "yellow")
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = "Clusters", ylab = "WGCNA Module Eigengene")
axis(1, at = 1:16, labels = 0:15)
matpoints(plotdat, col = modules, pch = 21)

cluster <- "turquoise"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, las = 3, at = 1:68, labels = unique(Idents(pbmc)))
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
abline(v=34.5,col="black")


tiff("/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/LglKD+CtrlWT_FC/WGCNA/WGCNA.tiff", height = 6000, width = 10000, res = 500)
par(mfrow = c(4,4))
cluster <- "1"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "2"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "3"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "4"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "5"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "6"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "7"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "8"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "9"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "10"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "11"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "12"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "13"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "14"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "15"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
cluster <- "16"
modules <- cluster
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = cluster, ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 0:17)
matpoints(plotdat, col = modules, pch = 21)
abline(h=0.000,col="black")
dev.off()

cluster <- colnames(datExpr)[net$colors == "1"]
write.csv(cluster, file = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/LglKD+CtrlWT_FC/WGCNA/1")
cluster <- colnames(datExpr)[net$colors == "2"]
write.csv(cluster, file = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/LglKD+CtrlWT_FC/WGCNA/2")
cluster <- colnames(datExpr)[net$colors == "3"]
write.csv(cluster, file = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/LglKD+CtrlWT_FC/WGCNA/3")
cluster <- colnames(datExpr)[net$colors == "4"]
write.csv(cluster, file = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/LglKD+CtrlWT_FC/WGCNA/4")
cluster <- colnames(datExpr)[net$colors == "5"]
write.csv(cluster, file = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/LglKD+CtrlWT_FC/WGCNA/5")
cluster <- colnames(datExpr)[net$colors == "6"]
write.csv(cluster, file = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/LglKD+CtrlWT_FC/WGCNA/6")
cluster <- colnames(datExpr)[net$colors == "7"]
write.csv(cluster, file = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/LglKD+CtrlWT_FC/WGCNA/7")
cluster <- colnames(datExpr)[net$colors == "8"]
write.csv(cluster, file = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/LglKD+CtrlWT_FC/WGCNA/8")
cluster <- colnames(datExpr)[net$colors == "9"]
write.csv(cluster, file = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/LglKD+CtrlWT_FC/WGCNA/9")



pbmc.data <- Read10X(data.dir = "/home/dchatterjee/Documents/Ctrl_disc/filtered/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "Ctrl_disc", min.cells = 3, min.features = 200)


##CCA
stim <- readRDS(file = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/NICD_FC/NICD_ONLYFC.rds")
ctrl <- readRDS(file = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/UpdLglIR_FC/UpdLglIR_LglIR_WT_onlyFC.rds")

stim@assays$spliced -> stim@assays$RNA
stim@meta.data$nCount_spliced -> stim@meta.data$nCount_RNA
stim@meta.data$nFeature_spliced -> stim@meta.data$nFeature_RNA
stim@meta.data$orig.ident -> stim@meta.data$sample
ctrl@meta.data$orig.ident <- ctrl@meta.data$sample

all.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:50)
all.combined <- IntegrateData(anchorset = all.anchors, dims = 1:50)
DefaultAssay(all.combined) <- "integrated"
all.combined <- ScaleData(all.combined, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)
all.combined <- RunPCA(all.combined, features = VariableFeatures(all.combined), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
all.combined <- FindNeighbors(all.combined, dims = 1:40)
all.combined <- FindClusters(all.combined, resolution = 0.5)
table(all.combined@active.ident,all.combined@meta.data$sample)
all.combined <- RunUMAP(all.combined, dims = 1:40, n.neighbors = 20, min.dist = 0.35)
DimPlot(all.combined, reduction = "umap", pt.size = 1, label = TRUE, group.by = "sample")
