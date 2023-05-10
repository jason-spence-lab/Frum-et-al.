library(Seurat)
library(SeuratWrappers)
library(CellChat)
library(slingshot)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
#Figure 1
#Integrate Fetal Data, Extract Epithelium, Identify Populations for Cell Chat Analysis
f59ddistal.data <- Read10X(data.dir = "./")
f59dairway.data <- Read10X(data.dir = "./")
f132ddistal.data <- Read10X(data.dir = "./")
f80ddistal.data <- Read10X(data.dir = "./")
f80dairway.data <- Read10X(data.dir = "./")
f103ddistal.data <- Read10X(data.dir = "./")
f103dairway.data <- Read10X(data.dir = "./")


f103dairway <- CreateSeuratObject(counts = f103dairway.data, project = "f103dairway", min.cells = 3, min.features = 200)
f103ddistal <- CreateSeuratObject(counts = f103ddistal.data, project = "f103ddistal", min.cells = 3, min.features = 200)
f132ddistal <- CreateSeuratObject(counts = f132ddistal.data, project = "f132ddistal", min.cells = 3, min.features = 200)
f59dairway <- CreateSeuratObject(counts = f59dairway.data, project = "f59dairway", min.cells = 3, min.features = 200)
f59ddistal <- CreateSeuratObject(counts = f59ddistal.data, project = "f59ddistal", min.cells = 3, min.features = 200)
f80dairway  <- CreateSeuratObject(counts = f80dairway.data, project = "f80dairway", min.cells = 3, min.features = 200)
f80ddistal <- CreateSeuratObject(counts = f80ddistal.data, project = "f80ddistal", min.cells = 3, min.features = 200)


f103dairway <- RenameCells(object = f103dairway, add.cell.id = "f103dairway")
f103ddistal <- RenameCells(object = f103ddistal, add.cell.id = "f103ddistal")
f132ddistal <- RenameCells(object = f132ddistal, add.cell.id = "f132ddistal")
f59dairway <- RenameCells(object = f59dairway, add.cell.id = "f59dairway")
f59ddistal <- RenameCells(object = f59ddistal, add.cell.id = "f59ddistal")
f80dairway  <- RenameCells(object = f80dairway, add.cell.id = "f80dairway")
f80ddistal <- RenameCells(object = f80ddistal, add.cell.id = "f80ddistal")


f103dairway[["percent.mt"]] <- PercentageFeatureSet(f103dairway, pattern = "^MT-") 
f103ddistal[["percent.mt"]] <- PercentageFeatureSet(f103ddistal, pattern = "^MT-") 
f132ddistal[["percent.mt"]] <- PercentageFeatureSet(f132ddistal, pattern = "^MT-") 
f59dairway[["percent.mt"]] <- PercentageFeatureSet(f59dairway, pattern = "^MT-") 
f59ddistal[["percent.mt"]] <- PercentageFeatureSet(f59ddistal, pattern = "^MT-") 
f80dairway[["percent.mt"]]  <- PercentageFeatureSet(f80dairway, pattern = "^MT-") 
f80ddistal[["percent.mt"]] <- PercentageFeatureSet(f80ddistal, pattern = "^MT-")

fetal.combined <- merge(f103dairway, y = c(f103ddistal, f132ddistal, f59dairway,f59ddistal, f80dairway, f80ddistal), project = "fetallung")

fetal.combined <- subset(fetal.combined, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

fetal.combined <- SplitObject(fetal.combined, split.by = "orig.ident")
for (i in seq_along(fetal.combined)) {
  fetal.combined[[i]] <- NormalizeData(fetal.combined[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(fetal.combined)
for (i in seq_along(along.with = fetal.combined)) {
  fetal.combined[[i]] <- ScaleData(fetal.combined[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(fetal.combined, anchor.features = features, reduction = "rpca", dims = 1:30)
fetal.combined.integrated <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(fetal.combined.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
fetal.combined.integrated <- CellCycleScoring(fetal.combined.integrated, s.features = s.genes, g2m.features = g2m.genes)
fetal.combined.integrated <- ScaleData(fetal.combined.integrated, vars.to.regress = c("S.Score", "G2M.Score"))
fetal.combined.integrated <- RunPCA(fetal.combined.integrated)
ElbowPlot(fetal.combined.integrated, ndims = 50)
fetal.combined.integrated <- RunUMAP(fetal.combined.integrated, dims = 1:18, reduction.name = "umap", return.model = TRUE)
fetal.combined.integrated <- FindNeighbors(fetal.combined.integrated, dims = 1:18)
fetal.combined.integrated <- FindClusters(fetal.combined.integrated, resolution = 0.5)
DimPlot(fetal.combined.integrated, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)

DefaultAssay(fetal.combined.integrated) <- "RNA"
fetal.combined.integrated <- NormalizeData(fetal.combined.integrated)

fetal.combined.simplified <- subset(fetal.combined.integrated, idents = c(0, 1, 5, 21, 19, 9, 7)) #subset populations of interest: Bud-Tip Progenitors (9), RSPO2+ Mesenchyme (0, 1), Differentiating Epithelium (7,21,19), Myofibroblasts (5)

DefaultAssay(fetal.combined.simplified) <- "integrated"
fetal.combined.simplified <- RunPCA(fetal.combined.simplified, verbose = FALSE)
ElbowPlot(fetal.combined.simplified, ndims = 50)
fetal.combined.simplified <- RunUMAP(fetal.combined.simplified, reduction = "pca", dims = 1:12)
fetal.combined.simplified <- FindNeighbors(fetal.combined.simplified, dims = 1:12)
fetal.combined.simplified <- FindClusters(fetal.combined.simplified, resolution = 0.15) 
DimPlot(fetal.combined.simplified, label = TRUE, pt.size = 2)

fetal.combined.simplified <- subset(fetal.combined.simplified, idents = c(0, 1, 2, 3)) #clean out weird small clusters
DefaultAssay(fetal.combined.simplified) <- "integrated"
fetal.combined.simplified <- RunPCA(fetal.combined.simplified, verbose = FALSE)
ElbowPlot(fetal.combined.simplified, ndims = 50)
fetal.combined.simplified <- RunUMAP(fetal.combined.simplified, reduction = "pca", dims = 1:10)
fetal.combined.simplified <- FindNeighbors(fetal.combined.simplified, dims = 1:10)
fetal.combined.simplified <- FindClusters(fetal.combined.simplified, resolution = 0.15) 
DimPlot(fetal.combined.simplified, label = TRUE, pt.size = 2)                                    

DefaultAssay(fetal.combined.simplified) <- "RNA"
fetal.combined.simplified <- NormalizeData(fetal.combined.simplified)   

fetal.combined.simplified.renamed <- RenameIdents(object = fetal.combined.simplified, '0' = "RSPO2+", '1' = "Myofibroblasts", '2' = "Differentiating Epithelium", '3' = "Bud Tip Progenitors")
fetal.combined.simplified.renamed.cellchat <- createCellChat(fetal.combined.simplified.renamed, group.by = "ident", assay = "RNA")
#CellChat Analysis
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling


fetal.combined.simplified.renamed.cellchat@DB <- CellChatDB.use
fetal.combined.simplified.renamed.cellchat <- subsetData(fetal.combined.simplified.renamed.cellchat)
fetal.combined.simplified.renamed.cellchat <- identifyOverExpressedGenes(fetal.combined.simplified.renamed.cellchat)
fetal.combined.simplified.renamed.cellchat <- identifyOverExpressedInteractions(fetal.combined.simplified.renamed.cellchat)
fetal.combined.simplified.renamed.cellchat <- projectData(fetal.combined.simplified.renamed.cellchat, PPI.human)
fetal.combined.simplified.renamed.cellchat <- computeCommunProb(fetal.combined.simplified.renamed.cellchat, raw.use = FALSE)


fetal.combined.simplified.renamed.cellchat <- filterCommunication(fetal.combined.simplified.renamed.cellchat, min.cells = 10)
fetal.combined.simplified.renamed.cellchat <- computeCommunProbPathway(fetal.combined.simplified.renamed.cellchat, thresh = 0.05)
fetal.combined.simplified.renamed.cellchat <- aggregateNet(fetal.combined.simplified.renamed.cellchat)
fetal.combined.simplified.renamed.cellchat <- netAnalysis_computeCentrality(fetal.combined.simplified.renamed.cellchat, slot.name = "netP")

groupSize <- as.numeric(table(fetal.combined.simplified.renamed.cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(fetal.combined.simplified.renamed.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(fetal.combined.simplified.renamed.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#visualize communication between populations
mat <- fetal.combined.simplified.renamed.cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
fetal.combined.simplified.renamed.cellchat@netP$pathways


fetal.combined.simplified.renamed.cellchat <- netAnalysis_computeCentrality(fetal.combined.simplified.renamed.cellchat, slot.name = "netP")

#Figure 1a
ppdf(file.path("./", paste0("Fetal TGFb", ".pdf")), w=6, h=6)
pathways.show <- c("TGFb")
netVisual_aggregate(fetal.combined.simplified.renamed.cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.label.cex = 0.001, arrow.size = 0.5, arrow.width = 2)
dev.off()
pdf(file.path("./", paste0("Fetal BMP", ".pdf")), w=6, h=6)
pathways.show <- c("BMP")
netVisual_aggregate(fetal.combined.simplified.renamed.cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.label.cex = 0.001, arrow.size = 0.5,arrow.width = 2)
dev.off()
pdf(file.path("./", paste0("Fetal Legend", ".pdf")), w=6, h=6)

netVisual_aggregate(fetal.combined.simplified.renamed.cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.label.cex = 0.001, arrow.size = 0.5,arrow.width = 2)
dev.off()
pdf(file.path("./", paste0("Fetal WNT", ".pdf")), w=6, h=6)
pathways.show <- c("WNT")
netVisual_aggregate(fetal.combined.simplified.renamed.cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,vertex.label.cex = 0.001, arrow.size = 0.5,arrow.width = 2)
dev.off()
pdf(file.path("./", paste0("Fetal FGF", ".pdf")), w=6, h=6)
pathways.show <- c("FGF")
netVisual_aggregate(fetal.combined.simplified.renamed.cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.label.cex = 0.001, arrow.size = 0.5,arrow.width = 2)
dev.off()

#Figure 1b 
#process each distal fetal dataset separately to identify and extract budtips

VlnPlot(f103ddistal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
f103ddistal <- subset(f103ddistal, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10)
s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes
f103ddistal <- CellCycleScoring(f103ddistal, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
f103ddistal <- SCTransform(f103ddistal, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
f103ddistal <- RunPCA(f103ddistal, features = VariableFeatures(object = f103ddistal))
ElbowPlot(f103ddistal, ndims = 50)
f103ddistal <- RunUMAP(f103ddistal, dims = 1:16, reduction.name = "umap", return.model = TRUE)
f103ddistal <- FindNeighbors(f103ddistal, dims = 1:16)
f103ddistal <- FindClusters(f103ddistal, resolution = 0.5)
DimPlot(f103ddistal, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)
pdf(file.path("./", paste0("FullData Compartment", ".pdf")), w=11, h=8.5)
FeaturePlot(f103ddistal, features = compartment, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f103ddistal, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f103ddistal, label = TRUE, pt.size = 0.5)
dev.off()
FeaturePlot(f103ddistal, features = compartment, pt.size = 0.2, order = TRUE)

f103ddistal.epithelial <- subset(f103ddistal, idents = c(6, 7))
DefaultAssay(f103ddistal.epithelial) <- "SCT"
f103ddistal.epithelial <- RunPCA(f103ddistal.epithelial, verbose = FALSE)
ElbowPlot(f103ddistal.epithelial, ndims = 50)
f103ddistal.epithelial <- RunUMAP(f103ddistal.epithelial, reduction = "pca", dims = 1:10)
f103ddistal.epithelial <- FindNeighbors(f103ddistal.epithelial, dims = 1:10)
f103ddistal.epithelial <- FindClusters(f103ddistal.epithelial, resolution = 0.5)
DimPlot(f103ddistal.epithelial, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)



pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f103ddistal.epithelial, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()


DefaultAssay(f103ddistal.epithelial) <- "RNA"
f103ddistal.epithelial <- NormalizeData(f103ddistal.epithelial)
pdf(file.path("./", paste0("Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f103ddistal.epithelial, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()

zero.markers <- FindMarkers(f103ddistal.epithelial, ident.1 = 0, min.pct = 0.25)
one.markers <- FindMarkers(f103ddistal.epithelial, ident.1 = 1, min.pct = 0.25)
two.markers <- FindMarkers(f103ddistal.epithelial, ident.1 = 2, min.pct = 0.25)
three.markers <- FindMarkers(f103ddistal.epithelial, ident.1 = 3, min.pct = 0.25)
four.markers <- FindMarkers(f103ddistal.epithelial, ident.1 = 4, min.pct = 0.25)

write.csv(zero.markers, "Cluster0.csv")
write.csv(one.markers, "Cluster1.csv")
write.csv(two.markers, "Cluster2.csv")
write.csv(three.markers, "Cluster3.csv")
write.csv(four.markers, "Cluster4.csv")

f103ddistal.budtip.cellids <- WhichCells(f103ddistal.epithelial, idents = c(1, 4, 2))
f103ddistal.bta.cellids <- WhichCells(f103ddistal.epithelial, idents = c(0))
f103ddistal.budtips <- subset(f103ddistal.epithelial, idents = c(1, 4, 2))
#f132distal
VlnPlot(f132ddistal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
f132ddistal <- subset(f132ddistal, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10)
s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes
f132ddistal <- CellCycleScoring(f132ddistal, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
f132ddistal <- SCTransform(f132ddistal, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
f132ddistal <- RunPCA(f132ddistal, features = VariableFeatures(object = f132ddistal))
ElbowPlot(f132ddistal, ndims = 50)
f132ddistal <- RunUMAP(f132ddistal, dims = 1:16, reduction.name = "umap", return.model = TRUE)
f132ddistal <- FindNeighbors(f132ddistal, dims = 1:16)
f132ddistal <- FindClusters(f132ddistal, resolution = 0.5)
DimPlot(f132ddistal, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)
pdf(file.path("./", paste0("FullData Compartment", ".pdf")), w=11, h=8.5)
FeaturePlot(f132ddistal, features = compartment, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f132ddistal, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f132ddistal, label = TRUE, pt.size = 0.5)
dev.off()
FeaturePlot(f132ddistal, features = compartment, pt.size = 0.2, order = TRUE)

f132ddistal.epithelial <- subset(f132ddistal, idents = c(7, 8))
DefaultAssay(f132ddistal.epithelial) <- "SCT"
f132ddistal.epithelial <- RunPCA(f132ddistal.epithelial, verbose = FALSE)
ElbowPlot(f132ddistal.epithelial, ndims = 50)
f132ddistal.epithelial <- RunUMAP(f132ddistal.epithelial, reduction = "pca", dims = 1:10)
f132ddistal.epithelial <- FindNeighbors(f132ddistal.epithelial, dims = 1:10)
f132ddistal.epithelial <- FindClusters(f132ddistal.epithelial, resolution = 0.5)

pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f132ddistal.epithelial, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()
pdf(file.path("./", paste0("UMAPbySample", ".pdf")), w=11, h=8.5)
DimPlot(f132ddistal.epithelial, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
dev.off()    

DefaultAssay(f132ddistal.epithelial) <- "RNA"
f132ddistal.epithelial <- NormalizeData(f132ddistal.epithelial)
pdf(file.path("./", paste0("Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f132ddistal.epithelial, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()


zero.markers <- FindMarkers(f132ddistal.epithelial, ident.1 = 0, min.pct = 0.25)
one.markers <- FindMarkers(f132ddistal.epithelial, ident.1 = 1, min.pct = 0.25)
two.markers <- FindMarkers(f132ddistal.epithelial, ident.1 = 2, min.pct = 0.25)
three.markers <- FindMarkers(f132ddistal.epithelial, ident.1 = 3, min.pct = 0.25)


write.csv(zero.markers, "Cluster0.csv")
write.csv(one.markers, "Cluster1.csv")
write.csv(two.markers, "Cluster2.csv")
write.csv(three.markers, "Cluster3.csv")



f132ddistal.budtip.cellids <- WhichCells(f132ddistal.epithelial, idents = c(1,3))
f132ddistal.bta.cellids <- WhichCells(f132ddistal.epithelial, idents = c(0))
f132ddistal.mixed.cellids <- WhichCells(f132ddistal.epithelial, idents = c(2))
f132ddistal.budtip <- subset(f132ddistal.epithelial, idents = c(1, 3))

#80 day distal fetal lung
VlnPlot(f80ddistal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
f80ddistal <- subset(f80ddistal, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10)
s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes
f80ddistal <- CellCycleScoring(f80ddistal, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
f80ddistal <- SCTransform(f80ddistal, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
f80ddistal <- RunPCA(f80ddistal, features = VariableFeatures(object = f80ddistal))
ElbowPlot(f80ddistal, ndims = 50)
f80ddistal <- RunUMAP(f80ddistal, dims = 1:22, reduction.name = "umap", return.model = TRUE)
f80ddistal <- FindNeighbors(f80ddistal, dims = 1:22)
f80ddistal <- FindClusters(f80ddistal, resolution = 0.5)
DimPlot(f80ddistal, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)
pdf(file.path("./", paste0("FullData Compartment", ".pdf")), w=11, h=8.5)
FeaturePlot(f80ddistal, features = compartment, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f80ddistal, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f80ddistal, label = TRUE, pt.size = 0.5)
dev.off()

f80ddistal.epithelial <- subset(f80ddistal, idents = c(7, 13, 10))
DefaultAssay(f80ddistal.epithelial) <- "SCT"
f80ddistal.epithelial <- RunPCA(f80ddistal.epithelial, verbose = FALSE)
ElbowPlot(f80ddistal.epithelial, ndims = 50)
f80ddistal.epithelial <- RunUMAP(f80ddistal.epithelial, reduction = "pca", dims = 1:10)
f80ddistal.epithelial <- FindNeighbors(f80ddistal.epithelial, dims = 1:10)
f80ddistal.epithelial <- FindClusters(f80ddistal.epithelial, resolution = 0.5)
DimPlot(f80ddistal.epithelial, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)
FeaturePlot(f80ddistal.epithelial, features = epithelial.markers, pt.size = 0.2, order = TRUE)

pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f80ddistal.epithelial, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()
pdf(file.path("./", paste0("UMAPbySample", ".pdf")), w=11, h=8.5)
DimPlot(f80ddistal.epithelial, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
dev.off()    

DefaultAssay(f80ddistal.epithelial) <- "RNA"
f80ddistal.epithelial <- NormalizeData(f80ddistal.epithelial)
pdf(file.path("./", paste0("Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f80ddistal.epithelial, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()




zero.markers <- FindMarkers(f80ddistal.epithelial, ident.1 = 0, min.pct = 0.25)
one.markers <- FindMarkers(f80ddistal.epithelial, ident.1 = 1, min.pct = 0.25)
two.markers <- FindMarkers(f80ddistal.epithelial, ident.1 = 2, min.pct = 0.25)
three.markers <- FindMarkers(f80ddistal.epithelial, ident.1 = 3, min.pct = 0.25)
four.markers <- FindMarkers(f80ddistal.epithelial, ident.1 = 4, min.pct = 0.25)


write.csv(zero.markers, "Cluster0.csv")
write.csv(one.markers, "Cluster1.csv")
write.csv(two.markers, "Cluster2.csv")
write.csv(three.markers, "Cluster3.csv")
write.csv(four.markers, "Cluster4.csv")

f80ddistal.budtip.cellids <- WhichCells(f80ddistal.epithelial, idents = c(0))
f80ddistal.bta.cellids <- WhichCells(f80ddistal.epithelial, idents = c(1))
f80ddistal.budtip <- subset(f80ddistal.epithelial, idents = c(0,3))

##59 day distal fetal lung
VlnPlot(f59ddistal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
f59ddistal <- subset(f59ddistal, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10)
s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes
f59ddistal <- CellCycleScoring(f59ddistal, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
f59ddistal <- SCTransform(f59ddistal, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
f59ddistal <- RunPCA(f59ddistal, features = VariableFeatures(object = f59ddistal))
ElbowPlot(f59ddistal, ndims = 50)
f59ddistal <- RunUMAP(f59ddistal, dims = 1:22, reduction.name = "umap", return.model = TRUE)
f59ddistal <- FindNeighbors(f59ddistal, dims = 1:22)
f59ddistal <- FindClusters(f59ddistal, resolution = 0.5)
DimPlot(f59ddistal, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)
pdf(file.path("./", paste0("FullData Compartment", ".pdf")), w=11, h=8.5)
FeaturePlot(f59ddistal, features = compartment, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f59ddistal, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f59ddistal, label = TRUE, pt.size = 0.5)
dev.off()

f59ddistal.epithelial <- subset(f59ddistal, idents = c(10, 13))
DefaultAssay(f59ddistal.epithelial) <- "SCT"
f59ddistal.epithelial <- RunPCA(f59ddistal.epithelial, verbose = FALSE)
ElbowPlot(f59ddistal.epithelial, ndims = 50)
f59ddistal.epithelial <- RunUMAP(f59ddistal.epithelial, reduction = "pca", dims = 1:10)
f59ddistal.epithelial <- FindNeighbors(f59ddistal.epithelial, dims = 1:10)
f59ddistal.epithelial <- FindClusters(f59ddistal.epithelial, resolution = 1.45)
DimPlot(f59ddistal.epithelial, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)

pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f59ddistal.epithelial, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()
pdf(file.path("./", paste0("UMAPbySample", ".pdf")), w=11, h=8.5)
DimPlot(f59ddistal.epithelial, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
dev.off()    

DefaultAssay(f59ddistal.epithelial) <- "RNA"
f59ddistal.epithelial <- NormalizeData(f59ddistal.epithelial)
pdf(file.path("./", paste0("Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f59ddistal.epithelial, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()


zero.markers <- FindMarkers(f59ddistal.epithelial, ident.1 = 0, min.pct = 0.25)
one.markers <- FindMarkers(f59ddistal.epithelial, ident.1 = 1, min.pct = 0.25)
two.markers <- FindMarkers(f59ddistal.epithelial, ident.1 = 2, min.pct = 0.25)
three.markers <- FindMarkers(f59ddistal.epithelial, ident.1 = 3, min.pct = 0.25)


write.csv(zero.markers, "Cluster0.csv")
write.csv(one.markers, "Cluster1.csv")
write.csv(two.markers, "Cluster2.csv")
write.csv(three.markers, "Cluster3.csv")


f59ddistal.budtip.cellids <- WhichCells(f59ddistal.epithelial, idents = c(1))
f59ddistal.bta.cellids <- WhichCells(f59ddistal.epithelial, idents = c(3))
f59ddistal.budtip <- subset(f59ddistal.epithelial, idents = c(1))

allbudtips <- merge(f59ddistal.budtip, c(f80ddistal.budtip, f103ddistal.budtips, f132ddistal.budtip))

allbudtips$orig.ident <- factor(allbudtips$orig.ident, levels = c("f132ddistal", "f103ddistal", "f80ddistal", "f59ddistal"))
DefaultAssay(allbudtips) <- "RNA"
allbudtips <- NormalizeData(allbudtips)


Dotplot_Zhiwei_Version <- function(seurat_object, gene_list) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 20, scale = TRUE, group.by = "orig.ident") + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8))
}

##Figure 1b
pdf(file.path("./", paste0("Fetal BudTip Only BMP4 ID2 SOX2 SOX9 dotplot and others", ".pdf")), w=8.0, h=7.5)
Dotplot_Zhiwei_Version(allbudtips, c("BMP2", "BMP4", "BMP5", "ID2" ,"TGFB2", "TGFB3", "SOX2", "SOX9", "SFTPC"))
dev.off()

pdf(file.path("./", paste0("Fetal BudTip Only BMP4 ID2 SOX2 SOX9 dotplot and others Aug 2022", ".pdf")), w=8.0, h=7.5)
Dotplot_Zhiwei_Version(allbudtips, c("BMP2", "BMP4", "BMP5","GDF5", "ID2" ,"TGFB1","TGFB2", "TGFB3", "SOX2", "SOX9", "SFTPC"))
dev.off()

##Figure 1c
rspo2 <- subset(fetal.combined.simplified.renamed, idents = "RSPO2+")
rspo2.59ddistal <- subset(rspo2, subset = orig.ident == "f59ddistal")
rspo2.80ddistal <- subset(rspo2, subset = orig.ident == "f80ddistal")
rspo2.103ddistal <- subset(rspo2, subset = orig.ident == "f103ddistal")
rspo2.132ddistal <- subset(rspo2, subset = orig.ident == "f132ddistal")

rspo2.merge <- merge(rspo2.59ddistal, c(rspo2.80ddistal, rspo2.103ddistal, rspo2.132ddistal))
DefaultAssay(rspo2.merge) <- "RNA"
rspo2.merge <- NormalizeData(rspo2.merge)

rspo2.merge$orig.ident <- factor(rspo2.merge$orig.ident, levels = c("f132ddistal", "f103ddistal", "f80ddistal", "f59ddistal"))

pdf(file.path("./", paste0("Fetal RSPO2 only BMP4 ID2 SOX2 SOX9 dotplot and others Aug 22", ".pdf")), w=8.0, h=7.5)
Dotplot_Zhiwei_Version(rspo2.merge, c("BMP2", "BMP4", "BMP5","GDF5", "ID2" ,"TGFB1", "TGFB2", "TGFB3", "RSPO2", "LIFR"))
dev.off()

#Explant Processing, Supplementary Figure 1f-g, Figure 1f-k. 

daythree.data <- Read10X(data.dir = "./")
daysix.data <- Read10X(data.dir = "./")
daynine.data <- Read10X(data.dir = "./")
daytwelve.data <- Read10X(data.dir = "./")


daythree <- CreateSeuratObject(counts = daythree.data, project = "daythree", min.cells = 3, min.features = 200)
daysix <- CreateSeuratObject(counts = daysix.data, project = "daysix", min.cells = 3, min.features = 200)
daynine <- CreateSeuratObject(counts = daynine.data, project = "daynine", min.cells = 3, min.features = 200)
daytwelve <- CreateSeuratObject(counts = daytwelve.data, project = "daytwelve", min.cells = 3, min.features = 200)

daythree <- RenameCells(object = daythree, add.cell.id = "daythree")
daysix <- RenameCells(object = daysix, add.cell.id = "daysix")
daynine <- RenameCells(object = daynine, add.cell.id = "daynine")
daytwelve <- RenameCells(object = daytwelve, add.cell.id = "daytwelve")


daythree[["percent.mt"]] <- PercentageFeatureSet(daythree, pattern = "^MT-") 
daysix[["percent.mt"]] <- PercentageFeatureSet(daysix, pattern = "^MT-") 
daynine[["percent.mt"]] <- PercentageFeatureSet(daynine, pattern = "^MT-") 
daytwelve[["percent.mt"]] <- PercentageFeatureSet(daytwelve, pattern = "^MT-")

explant.combined = merge(daythree, y = c(daysix, daynine, daytwelve))
explant.combined <- AddMetaData(explant.combined, "explant", col.name = "group")
VlnPlot(explant.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

explant.combined <- subset(explant.combined, subset = nFeature_RNA > 1000 & nFeature_RNA < 12000 & percent.mt < 10)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
explant.combined <- CellCycleScoring(explant.combined, s.features = s.genes, g2m.features = g2m.genes)
explant.combined.list <- SplitObject(explant.combined, split.by = "orig.ident")
explant.combined.list <- lapply(X = explant.combined.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
features <- SelectIntegrationFeatures(object.list = explant.combined.list)
explant.combined.list <- lapply(X = explant.combined.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = TRUE, vars.to.regress = c("S.Score", "G2M.Score"))
  x <- RunPCA(x, features = features, verbose = TRUE)
})
explant.anchors <- FindIntegrationAnchors(object.list = explant.combined.list, anchor.features = features, reduction = "rpca",  k.anchor = 5)
explant.combined.rpca <- IntegrateData(anchorset = explant.anchors)
DefaultAssay(explant.combined.rpca) <- "integrated"
explant.combined.rpca <- ScaleData(explant.combined.rpca, verbose = FALSE, vars.to.regress = c("S.Score", "G2M.Score"))
explant.combined.rpca <- RunPCA(explant.combined.rpca, npcs = 50, verbose = FALSE)
ElbowPlot(explant.combined.rpca, ndims = 50)
DefaultAssay(explant.combined.rpca) <- "integrated"

explant.combined.rpca <- RunUMAP(explant.combined.rpca, reduction = "pca", dims = 1:18)
explant.combined.rpca <- FindNeighbors(explant.combined.rpca, reduction = "pca", dims = 1:18)
explant.combined.rpca <- FindClusters(explant.combined.rpca, resolution = 0.5)

DimPlot(explant.combined.rpca, label = TRUE, pt.size = 2)

DefaultAssay(explant.combined.rpca) <- "RNA"
explant.combined.rpca <- NormalizeData(explant.combined.rpca)


explant.combined.rpca.renamed <- RenameIdents(explant.combined.rpca, '0' = "RSPO2+ Mesenchyme", '1' = "AT1", '2' = "Proliferative Epithelium", '3' = "Myofibroblast", '4' = "Pericyte", '5' = "Matrix Fibroblast", '6' = "Proliferative Fibroblast 1", '7' = "AT2", '8' = "Proliferative Fibroblast 2", '9' = "Endothelial", '10' = "Macrophages", '11' = "Airway", '12' = "Lymphatic", '13' = "Mesothelial", '14' = "T-Cell")

#Supplementary Figure 1f
pdf(file.path("./", paste0("Louvain no label", ".pdf")), w=11, h=8.5)
DimPlot(explant.combined.rpca, reduction = "umap", label = FALSE, pt.size = 2)
dev.off()

#Supplementary Figure 1g
levels(explant.combined.rpca) <- c("13","12", "9", "10", "14","8", "6", "5", "4","3", "0" ,"7", "2", "1", "11")
Dotplot_Zhiwei_Version <- function(seurat_object, gene_list) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 20, scale = TRUE) + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8))
}

pdf(file.path("./", paste0("Full Data Dotplot ", ".pdf")), w=18, h=11)
Dotplot_Zhiwei_Version(explant.combined.rpca, c("EPCAM", "CDH1", "SOX2", "AQP4",  "AGER", "SFTPC" , "LIFR","RSPO2", "PDGFRA", "ACTA2", "THY1", "PDGFRB","MGP", "TOP2A", "PCNA", "PTPRC", "CD3G","CD86", "PECAM1", "CA4","PROX1", "WT1", "UPK3B"))
dev.off()

#Data for Supplementary Figure 1h (graph made in prism)
write.csv(table(Idents(explant.combined.rpca), explant.combined.rpca$orig.ident))

#Supplementary Figure 1i
dayzero.data <- Read10X(data.dir = "./")
dayzero <- CreateSeuratObject(counts = dayzero.data, project = "dayzero", min.cells = 3, min.features = 200)
dayzero[["percent.mt"]] <- PercentageFeatureSet(dayzero, pattern = "^MT-")
setwd("~/scRNAseq/Output/2023-03-29 HT418 Explant Day 0/Filtered")
pdf(file.path("./", paste0("Sample QC", ".pdf")), w=11, h=8.5)
VlnPlot(dayzero, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
dayzero <- subset(dayzero, subset = nFeature_RNA > 1000 & nFeature_RNA < 12000 & percent.mt < 10)


dayzero <- SCTransform(dayzero)
DefaultAssay(dayzero) <- "SCT"
dayzero <- RunPCA(dayzero)
ElbowPlot(dayzero, ndims = 50)
dayzero <- RunUMAP(dayzero, dims = 1:20, reduction.name = "umap", return.model = TRUE)
dayzero <- FindNeighbors(dayzero, dims = 1:20)
dayzero <- FindClusters(dayzero, resolution = 0.5)
pdf(file.path("./", paste0("  Filtered Louvain", ".pdf")), w=11, h=8.5)
DimPlot(dayzero, label = TRUE, pt.size = 0.5)
dev.off()
DefaultAssay(dayzero) <- "RNA"
dayzero <- NormalizeData(dayzero)

dayzeroepithelium <- subset(dayzero, idents = c(5, 12))



setwd("~/scRNAseq/Output/2023-03-29 HT418 Explant Day 0/Epithelium Extracted/Filtered")



dayzeroepithelium <- SCTransform(dayzeroepithelium)
DefaultAssay(dayzeroepithelium) <- "SCT"
dayzeroepithelium <- RunPCA(dayzeroepithelium)
ElbowPlot(dayzeroepithelium, ndims = 50)
dayzeroepithelium <- RunUMAP(dayzeroepithelium, dims = 1:10, reduction.name = "umap", return.model = TRUE)
dayzeroepithelium <- FindNeighbors(dayzeroepithelium, dims = 1:10)
dayzeroepithelium <- FindClusters(dayzeroepithelium, resolution = 0.5)

dayzerobtp <- subset(dayzeroepithelium, idents = c(0, 3))
dayzerobtp <- AddMetaData(dayzerobtp, col.name = "class", "day 0 BTP")
explantclusterzero <- subset(explant.epithelium.rpca, idents = 0)
explantclusterzero <- AddMetaData(explantclusterzero, col.name = "class", "explant cluster 0")


explantclusterone <- subset(explant.epithelium.rpca, idents = 1)
explantclusterone <- AddMetaData(explantclusterone, col.name = "class", "explant cluster 1")




explantclustertwo <- subset(explant.epithelium.rpca, idents = 2)
explantclustertwo <- AddMetaData(explantclustertwo, col.name = "class", "explant cluster 2")

explantclusterthree <- subset(explant.epithelium.rpca, idents = 3)
explantclusterthree <- AddMetaData(explantclusterthree, col.name = "class", "explant cluster 3")

explantclusterfour <- subset(explant.epithelium.rpca, idents = 4)
explantclusterfour <- AddMetaData(explantclusterfour, col.name = "class", "explant cluster 4")

explantclusterfive <- subset(explant.epithelium.rpca, idents = 5)
explantclusterfive <- AddMetaData(explantclusterfive, col.name = "class", "explant cluster 5")

dayzerobtpvexplant <-merge(dayzerobtp, c(explantclusterzero, explantclusterone, explantclustertwo, explantclusterthree, explantclusterfour, explantclusterfive))
DefaultAssay(dayzerobtpvexplant) <- "RNA"
dayzerobtpvexplant <- NormalizeData(dayzerobtpvexplant)

dayzerobtpvexplant$class <- factor(dayzerobtpvexplant$class, levels = c("day 0 BTP", "explant cluster 0", "explant cluster 2", "explant cluster 3", "explant cluster 1", "explant cluster 4", "explant cluster 5"))


p <- VlnPlot(dayzerobtpvexplant, "SFTPC", group.by = "class", pt.size = 0, cols = c("#8494FF", "#F8766D", "#7CAE00", "#00BFC4", "#ABA300", "#00A9FF", "#FF61CC")) + NoLegend()
pdf(file.path("./", paste0("SFTPC Day 0 vs Clusters", ".pdf")), w=5.7, h=1.8)
p + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
dev.off()

p <- VlnPlot(dayzerobtpvexplant, "AGER", group.by = "class", pt.size = 0, cols = c("#8494FF", "#F8766D", "#7CAE00", "#00BFC4", "#ABA300", "#00A9FF", "#FF61CC")) + NoLegend()
pdf(file.path("./", paste0("AGER Day 0 vs Clusters", ".pdf")), w=5.7, h=1.8)
p + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
dev.off()


#Figure 1j
explant.epithelium.rpca <- subset(explant.combined.rpca, idents = c(1,7,2,11))

DefaultAssay(explant.epithelium.rpca) <- "integrated"
explant.epithelium.rpca <- RunPCA(explant.epithelium.rpca, verbose = FALSE)
ElbowPlot(explant.epithelium.rpca, ndims = 50)
explant.epithelium.rpca <- RunUMAP(explant.epithelium.rpca, reduction = "pca", dims = 1:7)
explant.epithelium.rpca <- FindNeighbors(explant.epithelium.rpca, dims = 1:7)
explant.epithelium.rpca <- FindClusters(explant.epithelium.rpca, resolution = 0.5)
DimPlot(explant.epithelium.rpca, label = TRUE, pt.size = 2)
epithelium.cluster5.cellids  <- WhichCells(epithelium.cluster5.cellids, idents = c(5)) ##remove cluster 5 which has low NKX2-1 expression

explant.epithelium.rpca <- subset(explant.epithelium.rpca, idents = c(0,1,2,3,4,6))
explant.epithelium.rpca <- RunPCA(explant.epithelium.rpca, verbose = FALSE)
ElbowPlot(explant.epithelium.rpca, ndims = 50)
explant.epithelium.rpca <- RunUMAP(explant.epithelium.rpca, reduction = "pca", dims = 1:7)
explant.epithelium.rpca <- FindNeighbors(explant.epithelium.rpca, dims = 1:7)
explant.epithelium.rpca <- FindClusters(explant.epithelium.rpca, resolution = 0.4)
pdf(file.path("./", paste0("Louvain no Label", ".pdf")), w=11, h=8.5)
DimPlot(explant.epithelium.rpca, reduction = "umap",  pt.size = 2)
dev.off()
#Supplementary Figure 1j
pdf(file.path("./", paste0("Louvain By Timepoint", ".pdf")), w=11, h=8.5)
DimPlot(explant.epithelium.rpca, reduction = "umap", group.by = "group", pt.size = 2)
dev.off()

#Data for Supplementary Figure 1k
write.csv(table(Idents(explant.epithelium.rpca), explant.epithelium.rpca$orig.ident))


#Figure 1k
DefaultAssay(explant.epithelium.rpca) <- "RNA"
explant.epithelium.rpca <- NormalizeData(explant.epithelium.rpca)

ExplantEpithelialCluster.markers <- c("SFTPC","SFTPB", "NAPSA", "ABCA3", "SLC34A2" ,"DMBT1","SFTPA1", "HOPX", "PDPN", "AGER", "RTKN2", "SPOCK2" ,"SCGB3A2", "SOX2","FOXJ1","TP63", "SCGB1A1", "SPDEF", "SOX9", "TESC", "SOX11", "CPM", "TOP2A", "MKI67")

explant.epithelium.rpca$seurat_clusters <- factor(explant.epithelium.rpca$seurat_clusters, levels = c("0", "5", "3", "2", "4", "1"))
pdf(file.path("./", paste0("Epithelial Marker DotPlot v2", ".pdf")), w=17, h=5)
Dotplot_Zhiwei_Version(explant.epithelium.rpca, ExplantEpithelialCluster.markers)
dev.off()


#Figure 1l
pdf(file.path("./", paste0("Ligands of interest dotplot.pdf")), w=6, h=5)
Dotplot_Zhiwei_Version(explant.epithelium.rpca, c("BMP2", "BMP4", "BMP5", "ID2" ,"TGFB2", "TGFB3")) 
dev.off()



#Supplementary Figure 1l

#pull AT2 cells from Travaglini et al., 2020 data
krasnow.data <- read.table(file=paste0("krasnow_hlca_10x_UMIs.csv"), sep = ",", row.names = 1, header = TRUE) #import data from Travaglini et al., 2020
krasnow <- CreateSeuratObject(krasnow.data, project = "krasnow", min.cells = 3, min.features = 200)
krasnow <- subset(krasnow, subset = nFeature_RNA > 500 & nFeature_RNA < 10000)
s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes
krasnow <- CellCycleScoring(krasnow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
krasnow <- SCTransform(krasnow, vars.to.regress = c("S.Score", "G2M.Score"))
krasnow <- RunPCA(krasnow, features = VariableFeatures(object = krasnow))
ElbowPlot(krasnow, ndims = 50)
krasnow <- RunUMAP(krasnow, dims = 1:25, reduction.name = "umap", return.model = TRUE)
krasnow <- FindNeighbors(krasnow, dims = 1:25)
krasnow <- FindClusters(krasnow, resolution = 0.5)

#extract epithelium from full data
krasnow.epithelial <- subset(krasnow, idents = c(16, 22, 14, 21, 4))
DefaultAssay(krasnow.epithelial) <- "SCT"
krasnow.epithelial <- RunPCA(krasnow.epithelial, verbose = FALSE)
ElbowPlot(krasnow.epithelial, ndims = 50)
krasnow.epithelial <- RunUMAP(krasnow.epithelial, reduction = "pca", dims = 1:20)
krasnow.epithelial <- FindNeighbors(krasnow.epithelial, dims = 1:20)
krasnow.epithelial <- FindClusters(krasnow.epithelial, resolution = 0.5)
DimPlot(krasnow.epithelial, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)

#extract at2 from epithelium
krasnow.at2 <- subset(krasnow.epithelial, idents = c(0, 1, 3))
krasnow.at2 <- AddMetaData(krasnow.at2, "Primary", col.name = "Sample")
#extract at2s from explant epithelium
explant.at2 <- subset(explant.epithelium.rpca, idents = c(1, 4))
explant.at2 <- AddMetaData(explant.at2, "Explant", col.name = "Sample")

#extract btps from fetal epithelium
fetal.btps <- merge(f59ddistal.budtip, f80ddistal.budtip, f103ddistal.budtips, f132distal.budtip)
fetal.btps <- AddMetaData(explant.at2, "Fetal", col.name = "Sample")
#merge adult AT2s with BTPs from fetal data and AT2s from explants
at2.combined <- merge(krasnow.at2, explant.at2, fetal.btps)

Dotplot_Zhiwei_Version <- function(seurat_object, gene_list) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 20, scale = TRUE, group.by = "Sample") + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8))
}


pdf(file.path("./", paste0("Fetal Explant Adult ATII marker DotPlot", ".pdf")), w=8, h=5)
Dotplot_Zhiwei_Version(at2.combined, c("SFTPC", "SFTPB",  "SLC34A2","MUC1", "NAPSA","LAMP3", "ABCA3", "SFTPA1"))
dev.off()

##Figure 3
frum.aec2.curated <- c("SFTPC", "SFTPB", "SFTPA1", "LAMP3", "NAPSA", "SLC34A2")
frum.btp.curated <- frum.btp.curated <- c("SOX9", "SOX2", "TESC", "SOX11", "CPM", "NKX2-1")

ckabday21.data <- Read10X(data.dir = "./")
threef.data <- Read10X(data.dir = "./") #BTP-Organoids under maintenance condition, 'day 0'
ckabday1.data <- Read10X(data.dir = "./")
ckabday6.data <- Read10X(data.dir = "./")


threef <- CreateSeuratObject(counts = threef.data, project = "threef", min.cells = 3, min.features = 200)
ckabday21 <- CreateSeuratObject(counts = ckabday21.data, project = "ckabday21", min.cells = 3, min.features = 200)
ckabday1 <- CreateSeuratObject(counts = ckabday1.data, project = "2fabday1", min.cells = 3, min.features = 200)
ckabday6 <- CreateSeuratObject(counts = ckabday6.data, project = "2fabday6", min.cells = 3, min.features = 200)
threef <- RenameCells(object = threef, add.cell.id = "3F")
threef[["percent.mt"]] <- PercentageFeatureSet(threef, pattern = "^MT-")
ckabday21 <- RenameCells(object = ckabday21, add.cell.id = "2FABDAY21")

ckabday1 <- RenameCells(object = ckabday1, add.cell.id = "2FABDAY1")
ckabday6 <- RenameCells(object = ckabday6, add.cell.id = "2FABDAY6")
ckabday1[["percent.mt"]] <- PercentageFeatureSet(ckabday1, pattern = "^MT-")
ckabday6[["percent.mt"]] <- PercentageFeatureSet(ckabday6, pattern = "^MT-")
ckabday21[["percent.mt"]] <- PercentageFeatureSet(ckabday21, pattern = "^MT-")


threef <- AddMetaData(threef, "Day 0", col.name = "group") 
ckabday6 <- AddMetaData(ckabday6, "Day 6", col.name = "group")
ckabday1 <- AddMetaData(ckabday1, "Day 1", col.name = "group")
ckabday21 <- AddMetaData(ckabday21, "Day 21", col.name = "group")

threef <- AddMetaData(threef, 0, col.name = "timepoint") 
ckabday6 <- AddMetaData(ckabday6, 1, col.name = "timepoint")
ckabday1 <- AddMetaData(ckabday1, 6, col.name = "timepoint")
ckabday21 <- AddMetaData(ckabday21, 21, col.name = "timepoint")


ckab <- merge(threef, y = c(ckabday1, ckabday6, ckabday21), project = "2FAB")

VlnPlot(ckab, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ckab <- subset(ckab, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)


DefaultAssay(ckab) <- "RNA"
ckab <- NormalizeData(ckab, verbose = TRUE)

#Figure 3c
ckab$group <- factor(ckab$group, levels = c("Day 0", "Day 1", "Day 6", "Day 21"))
p <- VlnPlot(ckab, features = frum.aec2.curated, group.by = "group",  ncol = 6, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[4]] <- p[[4]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[5]] <- p[[5]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[6]] <- p[[6]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN AT2 marker curated set", ".pdf")), w=5.7, h=1.8)
p
dev.off()

p <- VlnPlot(ckab, features = frum.btp.curated, group.by = "group",  ncol = 6, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[4]] <- p[[4]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[5]] <- p[[5]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[6]] <- p[[6]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN BTP marker curated set + NKX21", ".pdf")), w=5.7, h=1.8)
p
dev.off()


##Supplementary Figure 3
###Process each 2FAB timepoint individually and identify the best AT2s on the basis of AT2 marker expression
####Day 1
ckabday1 <- CreateSeuratObject(counts = ckabday1.data, project = "ckabday1", min.cells = 3, min.features = 200)
ckabday1 <- AddMetaData(ckabday1, "2FAB Day 1", col.name = "group")
ckabday1[["percent.mt"]] <- PercentageFeatureSet(ckabday1, pattern = "^MT-")
ckabday1 <- RenameCells(object = ckabday1, add.cell.id = "2FABDAY1")
VlnPlot(ckabday1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ckabday1 <- subset(ckabday1, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ckabday1 <- CellCycleScoring(ckabday1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ckabday1 <- SCTransform(ckabday1, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
ckabday1 <- RunPCA(ckabday1, features = VariableFeatures(object = ckabday1))

ElbowPlot(ckabday1, ndims = 50)
DefaultAssay(ckabday1) <- "SCT"
ckabday1 <- RunUMAP(ckabday1, dims = 1:18)
ckabday1 <- FindNeighbors(ckabday1, dims = 1:18)
ckabday1 <- FindClusters(ckabday1, resolution = 0.5)
DimPlot(ckabday1, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)


DefaultAssay(ckabday1) <- "RNA"
ckabday1 <- NormalizeData(ckabday1)

####Day 6
ckabday6 <- CreateSeuratObject(counts = ckabday6.data, project = "ckabday6", min.cells = 3, min.features = 200)
ckabday6 <- AddMetaData(ckabday6, "2FAB Day 6", col.name = "group")
ckabday6[["percent.mt"]] <- PercentageFeatureSet(ckabday6, pattern = "^MT-")
ckabday6 <- RenameCells(object = ckabday6, add.cell.id = "2FABDAY6")
VlnPlot(ckabday6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ckabday6 <- subset(ckabday6, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ckabday6 <- CellCycleScoring(ckabday6, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ckabday6 <- SCTransform(ckabday6, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
ckabday6 <- RunPCA(ckabday6, features = VariableFeatures(object = ckabday6))

ElbowPlot(ckabday6, ndims = 50)
DefaultAssay(ckabday6) <- "SCT"
ckabday6 <- RunUMAP(ckabday6, dims = 1:18)
ckabday6 <- FindNeighbors(ckabday6, dims = 1:18)
ckabday6 <- FindClusters(ckabday6, resolution = 0.5)
DimPlot(ckabday6, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)


DefaultAssay(ckabday6) <- "RNA"
ckabday6 <- NormalizeData(ckabday6)
####Day 21
ckabday21 <- CreateSeuratObject(counts = ckabday21.data, project = "ckabday21", min.cells = 3, min.features = 200)
ckabday21 <- AddMetaData(ckabday21, "2FAB Day 21", col.name = "group")
ckabday21[["percent.mt"]] <- PercentageFeatureSet(ckabday21, pattern = "^MT-")
ckabday21 <- RenameCells(object = ckabday21, add.cell.id = "2FABDAY21")
VlnPlot(ckabday21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ckabday21 <- subset(ckabday21, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ckabday21 <- CellCycleScoring(ckabday21, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ckabday21 <- SCTransform(ckabday21, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
ckabday21 <- RunPCA(ckabday21, features = VariableFeatures(object = ckabday21))

ElbowPlot(ckabday21, ndims = 50)
DefaultAssay(ckabday21) <- "SCT"
ckabday21 <- RunUMAP(ckabday21, dims = 1:18)
ckabday21 <- FindNeighbors(ckabday21, dims = 1:18)
ckabday21 <- FindClusters(ckabday21, resolution = 0.5)
DimPlot(ckabday21, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)

DefaultAssay(ckabday21) <- "RNA"
ckabday21 <- NormalizeData(ckabday21)

#Supplementary Figure 3e
pdf(file.path("./", paste0("Day 1 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabday1, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 1 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabday1, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 1 SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabday1, features =  "SFTPB" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 1 TOP2A", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabday1, features = "TOP2A", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

pdf(file.path("./", paste0("Day 1 Louvain", ".pdf")), w=11, h=8.5)
DimPlot(ckabday1, pt.size = 2)
dev.off()
pdf(file.path("./", paste0("Day 1 Louvain Best AT2s", ".pdf")), w=11, h=8.5)
DimPlot(ckabday1, cells.highlight = ckabday1.bestat2.cellids, pt.size = 2, sizes.highlight = 2, cols.highlight = "lightblue")
dev.off()


#Supplementary Figure 3f
pdf(file.path("./", paste0("Day 6 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabday6, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 6 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabday6, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 6 SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabday6, features =  "SFTPB" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 6 TOP2A", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabday6, features = "TOP2A", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

pdf(file.path("./", paste0("Day 6 Louvain", ".pdf")), w=11, h=8.5)
DimPlot(ckabday6, pt.size = 2)
dev.off()
pdf(file.path("./", paste0("Day 6 Louvain Best AT2s", ".pdf")), w=11, h=8.5)
DimPlot(ckabday6, cells.highlight = ckabday6.bestat2.cellids, pt.size = 2, sizes.highlight = 2, cols.highlight = "lightblue")
dev.off()

#Supplementary Figure 3g

pdf(file.path("./", paste0("Day 21 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabday21, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 21 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabday21, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 21 SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabday21, features =  "SFTPB" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 21 TOP2A", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabday21, features = "TOP2A", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 21 Louvain", ".pdf")), w=11, h=8.5)
DimPlot(ckabday21, pt.size = 2)
dev.off()
pdf(file.path("./", paste0("Day 21Louvain Best AT2s", ".pdf")), w=11, h=8.5)
DimPlot(ckabday21, cells.highlight = ckabday21.bestat2.cellids, pt.size = 2, sizes.highlight = 2, cols.highlight = "lightblue")
dev.off()

##Figure 3 Continued
###Integrate BTP-Organoids (Day 0) with Day 1, 6, and 21 
ckab.fastmnn <- ckab
set.seed(888)
DefaultAssay(ckab.fastmnn) <- "RNA"
ckab.fastmnn <- NormalizeData(ckab.fastmnn,
                                 normalization.method = "LogNormalize", scale.factor =10000)
ckab.fastmnn <-  FindVariableFeatures(ckab.fastmnn)
ckab.fastmnn <- RunFastMNN(object.list = SplitObject(ckab.fastmnn, split.by = "orig.ident"))
ElbowPlot(ckab.fastmnn, ndims = 50)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ckab.fastmnn <- CellCycleScoring(ckab.fastmnn, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ckab.fastmnn <- RunUMAP(ckab.fastmnn, dims = 1:30,reduction = "mnn")
ckab.fastmnn <- FindNeighbors(ckab.fastmnn, dims = 1:30,reduction = "mnn")
ckab.fastmnn <- FindClusters(ckab.fastmnn, resolution = 0.3, algorithm = 1)
ckab.fastmnn$group <- factor(ckab.fastmnn$group, levels = c("Day 0", "Day 1", "Day 6", "Day 21"))

DimPlot(ckab.fastmnn, group.by = c("group", "ident"), label = TRUE, pt.size = 2)

##PrctCellExpringGene and calc_helper function from : https://github.com/satijalab/seurat/issues/371#issuecomment-486384854
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}
##Data for Figure 3d (graph made in prism)
PrctCellExpringGene(ckab.fastmnn ,genes =c("SFTPC", "SFTPA1", "SFTPB"), group.by = "group")

##Figure 3e
pdf(file.path("./", paste0("2FAB FastMNN by Phase", ".pdf")), w=11, h=8.5)
DimPlot(ckab.fastmnn, reduction = "umap",  pt.size = 2, group.by = "Phase")
dev.off()
##Data for Figure 3f (graph made in prism)
write.csv(table(ckab.fastmnn$group, ckab.fastmnn$Phase), "RunbyPhase.csv")
##Figure 3g
pdf(file.path("./", paste0("2FAB FastMNN by Sample", ".pdf")), w=11, h=8.5)
DimPlot(ckab.fastmnn, reduction = "umap",  pt.size = 2, group.by = "group")
dev.off()
##Figure 3h
DefaultAssay(ckab.fastmnn) <- "RNA"
ckab.fastmnn <- NormalizeData(ckab.fastmnn)
pdf(file.path("./", paste0("2FAB FASTMNN SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(ckab.fastmnn, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("2FAB FASTMNN SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(ckab.fastmnn, features = "SFTPB", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("2FAB FASTMNN SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(ckab.fastmnn, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("2FAB FASTMNN LAMP3", ".pdf")), w=11, h=8.5)
FeaturePlot(ckab.fastmnn, features = "LAMP3", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("2FAB FASTMNN KRT5", ".pdf")), w=11, h=8.5)
FeaturePlot(ckab.fastmnn, features = "KRT5", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
##Figure 3i
pdf(file.path("./", paste0("2FAB FastMNN", ".pdf")), w=11, h=8.5)
DimPlot(ckab.fastmnn, reduction = "umap",  pt.size = 2, label = FALSE)
dev.off()
##Figure 3j
Dotplot_Zhiwei_Version <- function(seurat_object, gene_list) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 20, scale = TRUE) + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8))
}



levels(ckab.fastmnn) <- c("6", "1", "0", "5", "3" ,"2", "4")

ckab.fastmnn.markers <- c("SOX9", "SOX2", "TESC", "SFTPC","NAPSA", "SLC34A2","HOPX", "AGER", "PDPN", "TP63", "FOXJ1" , "SCGB1A1", "SCGB3A2" ,"SPDEF","MUC5AC", "MUC5B","ASCL1", "CHGA", "TOP2A", "PCNA")
pdf(file.path("./", paste0("2FAB FASTMNN Dotplot", ".pdf")), w=14, h=5.75)
Dotplot_Zhiwei_Version(ckab.fastmnn, ckab.fastmnn.markers)
dev.off()

##Figure 3k
ckab.fastmnn <- FindClusters(ckab.fastmnn, resolution = 1, algorithm = 1)

sce.ckab.nosct <- as.SingleCellExperiment(ckab.fastmnn)
sce.ckab.nosct <- slingshot(sce.ckab.nosct, clusterLabels = sce.ckab.nosct@colData@listData[["ident"]]
                               , reducedDim = "UMAP",
                               start.clus = "6", allow.breaks = FALSE, stretch = 0, omega = TRUE, extend = "n")
lnes.ckab <- getLineages(reducedDim(sce.ckab.nosct,"UMAP"), sce.ckab.nosct@colData@listData[["ident"]], start.clus = "6")
pt.ckab.nosct <- slingPseudotime(sce.ckab.nosct, na=TRUE)
pt.ckab.nosct[is.na(pt.ckab.nosct)] <- 0

colors <- pal[cut(pt.ckab.nosct[,2], breaks = 100)]

pdf(file.path("./", paste0("2FAB Psuedotime with Lineage", ".pdf")), w=11, h=8.5)

plot(reducedDim(sce.ckab.nosct, "UMAP"), col = colors, pch = 16, cex = 1)

plot(SlingshotDataSet(sce.ckab.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 1, add = TRUE)
plot(SlingshotDataSet(sce.ckab.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 3, add = TRUE)
plot(SlingshotDataSet(sce.ckab.nosct), lwd =6, type = 'lineage', col = c("gray"), linInd = 5, add = TRUE)
plot(SlingshotDataSet(sce.ckab.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 4, add = TRUE)
plot(SlingshotDataSet(sce.ckab.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 6, add = TRUE)
plot(SlingshotDataSet(sce.ckab.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 7, add = TRUE)
plot(SlingshotDataSet(sce.ckab.nosct), lwd=8, type = 'lineage', col = c("red"), linInd = 2, add = TRUE)

dev.off()
##Figure 3l
pdf(file.path("./", paste0("2FAB FastMNN All Highlight", ".pdf")), w=11, h=8.5)
DimPlot(ckab.fastmnn, cells.highlight = list(ckabday21.bestat2.cellids, ckabday6.bestat2.cellids, ckabday1.bestat2.cellids, threef.cellids), cols.highlight = c("purple", "blue", "green", "yellow"), pt.size = 2, sizes.highlight = 2, order = TRUE)
dev.off()
##Figure 3m
##Compare Primary AT2s (see Figure 5) to BTP organoids to identify genes expected to onset during alveolar differentiation

threef <- CreateSeuratObject(counts = threef.data, project = "threef", min.cells = 3, min.features = 200)
threef <- RenameCells(object = threef, add.cell.id = "3F")
threef[["percent.mt"]] <- PercentageFeatureSet(threef, pattern = "^MT-")
adultnof.data  <- Read10X(data.dir = "./") ##see figure 5
adultnof <- CreateSeuratObject(counts = threef.data, project = "adult", min.cells = 3, min.features = 200)
adultnof <- RenameCells(object = adultnof, add.cell.id = "adult")
adultnof[["percent.mt"]] <- PercentageFeatureSet(adultnof, pattern = "^MT-")
adultnof.bestat2 <- subset(adultnof, cells = adultnof.bestat2.cellids)

btptoat2.invitro <-merge(threef, y = adultnof.bestat2)
btptoat2.invitro <- subset(btptoat2.invitro, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 20)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
btptoat2.invitro <- CellCycleScoring(btptoat2.invitro, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
btptoat2.invitro <- SCTransform(btptoat2.invitro, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
setwd("~/scRNAseq/Output/2021-10-13 HT313 Long Term Alveolar Differentiation/2021-10-18 BTPtoAT2InVitro")
btptoat2.invitro <- RunPCA(btptoat2.invitro, features = VariableFeatures(object = btptoat2.invitro))
ElbowPlot(btptoat2.invitro, ndims = 50)
DefaultAssay(btptoat2.invitro) <- "SCT"

btptoat2.invitro <- RunUMAP(btptoat2.invitro, dims = 1:15)
btptoat2.invitro <- FindNeighbors(btptoat2.invitro, dims = 1:15)
btptoat2.invitro <- FindClusters(btptoat2.invitro, resolution = 0.1)
DimPlot(btptoat2.invitro, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)

DefaultAssay(btptoat2.invitro) <- "RNA"
btptoat2.invitro <- NormalizeData(btptoat2.invitro, verbose = FALSE)

Adultat2vBTP <- FindMarkers(btptoat2.invitro, ident.1 = 0, ident.2 = 1, min.pct = 0.25)
write.csv(Adultat2vBTP, "In vitro Diff.csv")
#top 200 primary AT2 organoids enirched from primary AT2 organoid v BTP organoid comparison
invitro.at2.genes.gained.nolist <- list(c("SFTPC",
                                          "SFTPA1",
                                          "SFTPA2",
                                          "SFTPB",
                                          "NAPSA",
                                          "SPINK5",
                                          "SLC34A2",
                                          "SERPIND1",
                                          "HPGD",
                                          "SCGB3A1",
                                          "PIGR",
                                          "AQP1",
                                          "LRRK2",
                                          "SFTA2",
                                          "HLA-DRA",
                                          "HHIP",
                                          "CEACAM6",
                                          "CD74",
                                          "AQP5",
                                          "VEPH1",
                                          "NUPR1",
                                          "FTL",
                                          "SFTPD",
                                          "SERPINA1",
                                          "MFSD2A",
                                          "HLA-B",
                                          "LAMP3",
                                          "MT-ND4L",
                                          "ADGRF5",
                                          "HLA-DPA1",
                                          "HLA-DPB1",
                                          "SLPI",
                                          "S100A9",
                                          "CD36",
                                          "SUSD2",
                                          "DBI",
                                          "XIST",
                                          "TMEM213",
                                          "FASN",
                                          "ALOX15B",
                                          "CTSH",
                                          "CCND2",
                                          "C3",
                                          "MALL",
                                          "ALPL",
                                          "RASGRF1",
                                          "SLC22A31",
                                          "S100A14",
                                          "HOPX",
                                          "SFTA1P",
                                          "PHLDA2",
                                          "CYB5A",
                                          "CTSD",
                                          "LMO3",
                                          "CD59",
                                          "C19orf33",
                                          "SLC39A8",
                                          "C2",
                                          "MT-CO2",
                                          "CACNA2D2",
                                          "LGI3",
                                          "SDR16C5",
                                          "HLA-DRB1",
                                          "MICAL2",
                                          "TFPI",
                                          "TNC",
                                          "S100A6",
                                          "POLR2L",
                                          "KRT7",
                                          "SNX30",
                                          "HIP1",
                                          "SFRP5",
                                          "ABCA3",
                                          "CAT",
                                          "MYO1B",
                                          "IL18",
                                          "MT-ND5",
                                          "MPZL2",
                                          "CHI3L1",
                                          "DUOX1",
                                          "CDC25B",
                                          "SDC4",
                                          "AP000357.2",
                                          "SFN",
                                          "RPS26",
                                          "SELENOW",
                                          "B2M",
                                          "SELENOP",
                                          "ARPC1B",
                                          "C1orf116",
                                          "FTH1",
                                          "CSTB",
                                          "MT-ND3",
                                          "LPCAT1",
                                          "DRAM1",
                                          "ETV1",
                                          "NFIC",
                                          "NNMT",
                                          "NQO1",
                                          "NPC2",
                                          "ITGB6",
                                          "RAB27A",
                                          "HSPH1",
                                          "EPHX1",
                                          "LGALS3",
                                          "MSMO1",
                                          "FBP1",
                                          "GPX4",
                                          "BCL2L1",
                                          "NFIX",
                                          "HLA-C",
                                          "BRI3",
                                          "ICAM1",
                                          "SCD",
                                          "HLA-DMA",
                                          "C16orf89",
                                          "FGGY",
                                          "AK1",
                                          "FOLR1",
                                          "MSN",
                                          "ISCU",
                                          "HPCAL1",
                                          "LIPH",
                                          "MGST1",
                                          "CAPN8",
                                          "IL1R1",
                                          "PLXND1",
                                          "MEGF6",
                                          "MBNL1",
                                          "GGT5",
                                          "PDXK",
                                          "CXCL17",
                                          "CA2",
                                          "MT-ND4",
                                          "HIF1A",
                                          "TMEM238",
                                          "UQCR10",
                                          "MT-ATP6",
                                          "SCNN1A",
                                          "HHIP-AS1",
                                          "STC1",
                                          "CBR1",
                                          "MTUS1",
                                          "ACSL4",
                                          "RAB27B",
                                          "SNHG7",
                                          "LY6E",
                                          "MEGF9",
                                          "CEBPA",
                                          "SLC6A20",
                                          "SNX25",
                                          "ANXA1",
                                          "CYP1B1",
                                          "SDC1",
                                          "TGFBR2",
                                          "ADIPOR1",
                                          "DCXR",
                                          "CAPN2",
                                          "MMP28",
                                          "TSTD1",
                                          "ACSS2",
                                          "MT-CO1",
                                          "GGTLC1",
                                          "DCBLD2",
                                          "TGM2",
                                          "HSPB8",
                                          "CPM",
                                          "PID1",
                                          "RHOBTB2",
                                          "EPDR1",
                                          "QPRT",
                                          "SULT1A1",
                                          "ADGRF1",
                                          "TAOK3",
                                          "MID1IP1",
                                          "AQP4",
                                          "ALDH2",
                                          "DGKD",
                                          "ZMAT3",
                                          "CD9",
                                          "GALNT10",
                                          "MEG3",
                                          "IFI16",
                                          "CREG1",
                                          "HLA-A",
                                          "SELENBP1",
                                          "DHCR24",
                                          "ARHGDIB",
                                          "TST",
                                          "STOM",
                                          "ACSL1",
                                          "SLC66A1L",
                                          "CTSS",
                                          "FAH",
                                          "PARM1",
                                          "COX17",
                                          "TMEM125",
                                          "HLA-DOA",
                                          "MCUR1"))

#gene set module scoring of BTP-organoids and all CK + AB treatment timepoints
ckab <- merge(threef, y = c(ckabday1, ckabday6, ckabday21), project = "2FAB")
VlnPlot(ckab, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ckab <- subset(ckab, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
DefaultAssay(ckab) <- "RNA"
ckab <- NormalizeData(ckab, verbose = TRUE)
y <- ckab
y$group <- factor(y$group, levels = c("Day 0", "Day 1", "Day 6", "Day 21"))
y <- AddModuleScore(y, features = invitro.at2.genes.gained.nolist, name = "At2invitro", search = TRUE, assay = "RNA")
pdf(file.path("./", paste0("In Vitro Alveolar Differentation Scoring")), w=11, h=8.5)
print(VlnPlot(y, features = "At2invitro1", pt.size = 0, group.by = "group"))
dev.off()

##Figure 5
ckabnof.data <- Read10X(data.dir = "./") #CK + AB induced organoids for 21 days, expanded in  SFFFwithout FGF10 for 120 days

ckabnof <- CreateSeuratObject(counts = ckabnof.data, project = "ckabnof", min.cells = 3, min.features = 200)
ckabnof <- RenameCells(object = ckabnof, add.cell.id = "ckabnof")
ckabnof[["percent.mt"]] <- PercentageFeatureSet(ckabnof, pattern = "^MT-")

VlnPlot(ckabnof, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ckabnof <- subset(ckabnof, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 20)


ckabnof <- SCTransform(ckabnof)
ckabnof <- RunPCA(ckabnof)
ElbowPlot(ckabnof, ndims = 50)
ckabnof <- RunUMAP(ckabnof, dims = 1:16, reduction.name = "umap")
ckabnof <- FindNeighbors(ckabnof, dims = 1:16)
ckabnof <- FindClusters(ckabnof, resolution = 0.5)
DimPlot(ckabnof, label = TRUE)

##Figure 5a
pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(ckabnof, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(ckabnof) <- "RNA"
ckabnof <- NormalizeData(ckabnof)
ckabnof  <- CellCycleScoring(ckabnof, s.features = s.genes, g2m.features = g2m.genes)

PrctCellExpringGene(ckabnof, "MUC5AC")
FeaturePlot(ckabnof, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))

pdf(file.path("./", paste0("Day 120 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabnof, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 120 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabnof, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 120 SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabnof, features =  "SFTPB" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 120 TOP2A", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabnof, features = "TOP2A", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 120 MUC5AC", ".pdf")), w=11, h=8.5)
FeaturePlot(ckabnof, features = "MUC5AC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

#Figure 5b
pdf(file.path("./", paste0("AT2 Marker Marker Dotplot", ".pdf")), w=11, h=8.5)
Dotplot_Zhiwei_Version(ckabnof, c("SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "PGC", "NAPSA", "SFTA2", "SFTA3", "MUC5AC"))
dev.off()

#Supplementary Figure 5
adultnof.data  <- Read10X(data.dir = "./") #primary AT2 organoids cultured in  SFFFfor 90 days and then switched to primary AT2 organoid with no FGF10 for 30 days before sequencing

adultnof <- CreateSeuratObject(counts = adultnof.data, project = "adultnof", min.cells = 3, min.features = 200)
adultnof <- RenameCells(object = adultnof, add.cell.id = "adultnof")
adultnof[["percent.mt"]] <- PercentageFeatureSet(adultnof, pattern = "^MT-")

VlnPlot(adultnof, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
adultnof <- subset(adultnof, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 20)


adultnof <- SCTransform(adultnof)
adultnof <- RunPCA(adultnof)
ElbowPlot(adultnof, ndims = 50)
adultnof <- RunUMAP(adultnof, dims = 1:16, reduction.name = "umap")
adultnof <- FindNeighbors(adultnof, dims = 1:16)
adultnof <- FindClusters(adultnof, resolution = 0.5)
DimPlot(adultnof, label = TRUE)

##Supplementary Figure 5a
pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(adultnof, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(adultnof) <- "RNA"
adultnof <- NormalizeData(adultnof)



pdf(file.path("./", paste0("Adult NOF SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(adultnof, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult NOF SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(adultnof, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult NOF SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(adultnof, features =  "SFTPB" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult NOF TOP2A", ".pdf")), w=11, h=8.5)
FeaturePlot(adultnof, features = "TOP2A", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult NOF MUC5AC", ".pdf")), w=11, h=8.5)
FeaturePlot(adultnof, features = "MUC5AC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

#Supplementary Figure 5b
pdf(file.path("./", paste0("AT2 Marker Marker Dotplot", ".pdf")), w=11, h=8.5)
Dotplot_Zhiwei_Version(adultnof, c("SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "PGC", "NAPSA", "SFTA2", "SFTA3", "MUC5AC"))
dev.off()

##Supplementary Figure 5c
ckdcinof.data  <- Read10X(data.dir = "./") #CKDCI induced organoids for 21 days, expanded in  SFFFwithout FGF10 for 120 days
ckdcinof <- CreateSeuratObject(counts = ckdcinof.data, project = "ckdcinof", min.cells = 3, min.features = 200)
ckdcinof <- RenameCells(object = ckdcinof, add.cell.id = "ckdcinof")
ckdcinof[["percent.mt"]] <- PercentageFeatureSet(ckdcinof, pattern = "^MT-")

VlnPlot(ckdcinof, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ckdcinof <- subset(ckdcinof, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 20)
setwd("~/scRNAseq/Output/2022-03-16 Frum 2022 Figure 5/ckdci Nof")


ckdcinof <- SCTransform(ckdcinof)
ckdcinof <- RunPCA(ckdcinof)
ElbowPlot(ckdcinof, ndims = 50)
ckdcinof <- RunUMAP(ckdcinof, dims = 1:16, reduction.name = "umap")
ckdcinof <- FindNeighbors(ckdcinof, dims = 1:16)
ckdcinof <- FindClusters(ckdcinof, resolution = 0.5)
DimPlot(ckdcinof, label = TRUE)
write.csv(table(Idents(ckdcinof), ckdcinof$orig.ident), "ClusterbyRun.csv") #number of cells in each cluster by run

pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(ckdcinof, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(ckdcinof) <- "RNA"
ckdcinof <- NormalizeData(ckdcinof)

FeaturePlot(ckdcinof, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))

pdf(file.path("./", paste0("CKDCI NOF SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdcinof, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("CKDCI NOF SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdcinof, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("CKDCI NOF SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdcinof, features =  "SFTPB" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("CKDCI NOF TOP2A", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdcinof, features = "TOP2A", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("CKDCI NOF MUC5AC", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdcinof, features = "MUC5AC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

##Supplementary Figure 5d
pdf(file.path("./", paste0("AT2 Marker Marker Dotplot", ".pdf")), w=11, h=8.5)
Dotplot_Zhiwei_Version(ckdcinof, c("SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "PGC", "NAPSA", "SFTA2", "SFTA3", "MUC5AC"))
dev.off()

#Supplementary Figure 5f
day120.merge <- merge(ckabnof, c(ckdcinof, adultnof))
DefaultAssay(day120.merge) <- "RNA"
day120.merge <- NormalizeData(day120.merge)

otherlineages.markers <- c("AGER", "TP63", "FOXJ1" , "SCGB1A1", "ASCL1", "CHGA", "SPDEF", "MUC5AC")

VlnPlot(day120.merge, otherlineages.markers, group.by = "orig.ident")
p <- VlnPlot(day120.merge, features = otherlineages.markers, group.by = "orig.ident",  ncol = 4, pt.size = 0.5)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[4]] <- p[[4]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[5]] <- p[[5]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[6]] <- p[[6]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[7]] <- p[[7]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[8]] <- p[[8]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())

pdf(file.path("./", paste0("VLN Bulk Day 120 Alternative CellTypes Pt", ".pdf")), w=5.7, h=3.5)
p
dev.off()

p <- VlnPlot(day120.merge, features = otherlineages.markers, group.by = "orig.ident",  ncol = 4, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[4]] <- p[[4]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[5]] <- p[[5]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[6]] <- p[[6]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[7]] <- p[[7]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[8]] <- p[[8]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN Bulk Day 120 Alternative CellTypes", ".pdf")), w=5.7, h=3.5)
p
dev.off()

##Data for Figure 5c (graph made in prism)
PrctCellExpringGene(day120.merge, "SFTPC", group.by = "orig.ident")
PrctCellExpringGene(day120.merge, "SFTPA1", group.by = "orig.ident")
PrctCellExpringGene(day120.merge, "SFTPB", group.by = "orig.ident")

##Figure 5d
day120.merge$orig.ident <- factor(day120.merge$orig.ident, levels = c("adultnof", "ckabnof", "ckdcinof"))

VlnPlot(day120.merge, c("SFTPC", "SFTPA1", "SFTPB"), pt.size = 0, group.by = "orig.ident")

p <- VlnPlot(day120.merge, features = c("SFTPC", "SFTPA1", "SFTPB"), group.by = "orig.ident",  ncol = 4, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())

pdf(file.path("./", paste0("VLN Bulk Compare SFTPC A1 B", ".pdf")), w=5.7, h=1.5)
p
dev.off()

##Figure 5e
PrctCellExpringGene(day120.merge, "MUC5AC", group.by = "orig.ident")

##Figure 5f - h
##Reference Based Mapping of organoid samples to single cell sequencing of microdissected distal and proximal lung (Murthy et al., 2022)
#load Murthy et al., 2022
epi.sub_SAE10x3.integrated_v2 <- readRDS("~/filepath/epi.sub_SAE10x3.integrated_v2.RDS")
epi.sub_SAE10x3.integrated_v2 <- UpdateSeuratObject(epi.sub_SAE10x3.integrated_v2)
epi.sub_SAE10x3.integrated_v2 <- RunUMAP(object = epi.sub_SAE10x3.integrated_v2, dims = 1:20,   min.dist = 1, n.neighbors = 20, spread = 1.5, local.connectivity = 10, return.model = TRUE)
DefaultAssay(epi.sub_SAE10x3.integrated_v2) <- "RNA" #data already normalized

#Map primary AT2 organoids cultured in  SFFFfor 90 days and then switched to primary AT2 organoid with no FGF10 for 30 days before sequencing
adultnof.data  <- Read10X(data.dir = "./")
adultnof <- CreateSeuratObject(counts = adultnof.data, project = "adultnof", min.cells = 3, min.features = 200)
adultnof <- RenameCells(object = adultnof, add.cell.id = "adultnof")
adultnof[["percent.mt"]] <- PercentageFeatureSet(adultnof, pattern = "^MT-")
VlnPlot(adultnof, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
adultnof <- subset(adultnof, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 20)
adultnof <- NormalizeData(adultnof)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
adultnof <- CellCycleScoring(adultnof, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
adultnof <- FindVariableFeatures(adultnof, selection.method = "vst", nfeatures = 2000)
adultnof <- ScaleData(adultnof)
adultnof <- RunPCA(adultnof, features = VariableFeatures(adultnof))
ElbowPlot(adultnof, ndims = 50)
adultnof <- RunUMAP(adultnof, dims = 1:20, reduction.name = "umap", return.model = TRUE)
adultnof <- FindNeighbors(adultnof, dims = 1:20)
adultnof <- FindClusters(adultnof, resolution = 0.3)
DimPlot(adultnof, label = TRUE, pt.size = 1)

common.features <- intersect(rownames(epi.sub_SAE10x3.integrated_v2), rownames(adultnof))
length(x = common.features)
epi.sub_SAE10x3.integrated_v2.anchors <- FindTransferAnchors(
  reference = epi.sub_SAE10x3.integrated_v2,
  normalization.method = "LogNormalize",
  query = adultnof,
  k.filter = NA,
  reference.reduction = "pca",
  features = common.features,
  dims = 1:20
)


adultnof.mapped.tata <- MapQuery(
  anchorset = epi.sub_SAE10x3.integrated_v2.anchors,
  query = adultnof,
  reference = epi.sub_SAE10x3.integrated_v2,
  refdata = "cell_type",
  reference.reduction = "pca", 
  transferdata.args = list(k.weight = 10),
  integrateembeddings.args = list(k.weight = 10),
  reduction.model = "umap"
)
adultnof.label.at2 <- subset(adultnof.mapped.tata, predicted.id == "AT2")
adultnof.label.at2.cellids <- WhichCells(adultnof.label.at2)
adultnof.label.pat2 <- subset(adultnof.mapped.tata, predicted.id == "Proliferating AT2")
adultnof.label.pat2.cellids <- WhichCells(adultnof.label.pat2)
adultnof.label.ne <- subset(adultnof.mapped.tata, predicted.id == "Neuroendocrine")
adultnof.label.ne.cellids <- WhichCells(adultnof.label.ne)
adultnof.label.1a1neg <- subset(adultnof.mapped.tata, predicted.id == "SFTPB+ SCGB3A2+ SCGB1A1-")
adultnof.label.1a1neg.cellids <- WhichCells(adultnof.label.1a1neg)
adultnof.label.trb <- subset(adultnof.mapped.tata, predicted.id == "SFTPC+ SCGB3A2+")
adultnof.label.trb.cellids <- WhichCells(adultnof.label.trb)

adultnof.mapped.tata@reductions[["umap"]] <- adultnof.mapped.tata@reductions[["ref.umap"]]
adultnof.visual <- merge(epi.sub_SAE10x3.integrated_v2, adultnof.mapped.tata, merge.dr = "umap")

predictions.adultnof.tata <- table(adultnof.mapped.tata$predicted.id)

#Map CK + AB induced organoids for 21 days, expanded in  SFFFwithout FGF10 for 120 days
ckabnof.data <- Read10X(data.dir = "./")
ckabnof <- CreateSeuratObject(counts = ckabnof.data, project = "ckabnof", min.cells = 3, min.features = 200)
ckabnof <- RenameCells(object = ckabnof, add.cell.id = "ckabnof")
ckabnof[["percent.mt"]] <- PercentageFeatureSet(ckabnof, pattern = "^MT-")
VlnPlot(ckabnof, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ckabnof <- subset(ckabnof, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 20)

setwd("~/scRNAseq/Output/2022-03-16 Frum 2022 Figure 5/Label Transfer/ckabnof")

ckabnof <- NormalizeData(ckabnof)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ckabnof <- CellCycleScoring(ckabnof, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ckabnof <- FindVariableFeatures(ckabnof, selection.method = "vst", nfeatures = 2000)
ckabnof <- ScaleData(ckabnof)
ckabnof <- RunPCA(ckabnof, features = VariableFeatures(ckabnof))
ElbowPlot(ckabnof, ndims = 50)
ckabnof <- RunUMAP(ckabnof, dims = 1:20, reduction.name = "umap", return.model = TRUE)
ckabnof <- FindNeighbors(ckabnof, dims = 1:20)
ckabnof <- FindClusters(ckabnof, resolution = 0.3)
DimPlot(ckabnof, label = TRUE, pt.size = 1)

common.features <- intersect(rownames(epi.sub_SAE10x3.integrated_v2), rownames(ckabnof))
length(x = common.features)
epi.sub_SAE10x3.integrated_v2.anchors <- FindTransferAnchors(
  reference = epi.sub_SAE10x3.integrated_v2,
  normalization.method = "LogNormalize",
  query = ckabnof,
  k.filter = NA,
  reference.reduction = "pca",
  features = common.features,
  dims = 1:20
)


ckabnof.mapped.tata <- MapQuery(
  anchorset = epi.sub_SAE10x3.integrated_v2.anchors,
  query = ckabnof,
  reference = epi.sub_SAE10x3.integrated_v2,
  refdata = "cell_type",
  reference.reduction = "pca", 
  transferdata.args = list(k.weight = 10),
  integrateembeddings.args = list(k.weight = 10),
  reduction.model = "umap"
)

ckabnof.label.at2 <- subset(ckabnof.mapped.tata, predicted.id == "AT2")
ckabnof.label.at2.cellids <- WhichCells(ckabnof.label.at2)
ckabnof.label.at1 <- subset(ckabnof.mapped.tata, predicted.id == "AT1")
ckabnof.label.at1.cellids <- WhichCells(ckabnof.label.at1)
ckabnof.label.pat2 <- subset(ckabnof.mapped.tata, predicted.id == "Proliferating AT2")
ckabnof.label.pat2.cellids <- WhichCells(ckabnof.label.pat2)
ckabnof.label.diffbas <- subset(ckabnof.mapped.tata, predicted.id == "Differentiating Basal")
ckabnof.label.diffbas.cellids <- WhichCells(ckabnof.label.diffbas)
ckabnof.label.ne <- subset(ckabnof.mapped.tata, predicted.id == "Neuroendocrine")
ckabnof.label.ne.cellids <- WhichCells(ckabnof.label.ne)
ckabnof.label.1a1neg <- subset(ckabnof.mapped.tata, predicted.id == "SFTPB+ SCGB3A2+ SCGB1A1-")
ckabnof.label.1a1neg.cellids <- WhichCells(ckabnof.label.1a1neg)
ckabnof.label.1a1pos <- subset(ckabnof.mapped.tata, predicted.id == "SFTPB+ SCGB3A2+ SCGB1A1+")
ckabnof.label.1a1pos.cellids <- WhichCells(ckabnof.label.1a1pos)
ckabnof.label.goblet <- subset(ckabnof.mapped.tata, predicted.id == "MUC5AC+ MUC5B+")
ckabnof.label.goblet.cellids <- WhichCells(ckabnof.label.goblet)
ckabnof.mapped.tata@reductions[["umap"]] <- ckabnof.mapped.tata@reductions[["ref.umap"]]
ckabnof.visual <- merge(epi.sub_SAE10x3.integrated_v2, ckabnof.mapped.tata, merge.dr = "umap")



predictions.ckabnof.tata <- table(ckabnof.mapped.tata$predicted.id)


#Map CKDCI induced organoids for 21 days, expanded in  SFFFwithout FGF10 for 120 days
ckdcinof.data  <- Read10X(data.dir = "./")
ckdcinof <- CreateSeuratObject(counts = ckdcinof.data, project = "ckdcinof", min.cells = 3, min.features = 200)
ckdcinof <- RenameCells(object = ckdcinof, add.cell.id = "ckdcinof")
ckdcinof[["percent.mt"]] <- PercentageFeatureSet(ckdcinof, pattern = "^MT-")
VlnPlot(ckdcinof, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ckdcinof <- subset(ckdcinof, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 20)



ckdcinof <- NormalizeData(ckdcinof)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ckdcinof <- CellCycleScoring(ckdcinof, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ckdcinof <- FindVariableFeatures(ckdcinof, selection.method = "vst", nfeatures = 2000)
ckdcinof <- ScaleData(ckdcinof)
ckdcinof <- RunPCA(ckdcinof, features = VariableFeatures(ckdcinof))
ElbowPlot(ckdcinof, ndims = 50)
ckdcinof <- RunUMAP(ckdcinof, dims = 1:20, reduction.name = "umap", return.model = TRUE)
ckdcinof <- FindNeighbors(ckdcinof, dims = 1:20)
ckdcinof <- FindClusters(ckdcinof, resolution = 0.3)
DimPlot(ckdcinof, label = TRUE, pt.size = 1)

common.features <- intersect(rownames(epi.sub_SAE10x3.integrated_v2), rownames(ckdcinof))
length(x = common.features)
epi.sub_SAE10x3.integrated_v2.anchors <- FindTransferAnchors(
  reference = epi.sub_SAE10x3.integrated_v2,
  normalization.method = "LogNormalize",
  query = ckdcinof,
  k.filter = NA,
  reference.reduction = "pca",
  features = common.features,
  dims = 1:20
)


ckdcinof.mapped.tata <- MapQuery(
  anchorset = epi.sub_SAE10x3.integrated_v2.anchors,
  query = ckdcinof,
  reference = epi.sub_SAE10x3.integrated_v2,
  refdata = "cell_type",
  reference.reduction = "pca", 
  transferdata.args = list(k.weight = 10),
  integrateembeddings.args = list(k.weight = 10),
  reduction.model = "umap"
)




ckdcinof.mapped.tata@reductions[["umap"]] <- ckdcinof.mapped.tata@reductions[["ref.umap"]]
ckdcinof.visual <- merge(epi.sub_SAE10x3.integrated_v2, ckdcinof.mapped.tata, merge.dr = "umap")

predictions.ckdcinof.tata <- table(ckdcinof.mapped.tata$predicted.id)

##Figure 5F
cols <- c("AT1" = "#f5fc86", "AT2" = "#7CAE00", "Proliferating" = "#ED68ED", "Proliferating AT2" ="#FF61CC", "Differentiating Basal" = "#8494FF", "Neuroendocrine" = "#ab467a" , "SFTPB+ SCGB3A2+ SCGB1A1-" = "#ABA300", "SFTPB+ SCGB3A2+ SCGB1A1+" = "#E68613", "SFTPB- KRT5+ Basal" = "#00C19A",  "Proliferating" = "#ED68ED", "MUC5AC+ MUC5B+" = "#00BE67", "Immature AT1" ="#0CB702", "Ciliated" ="#F8766D", "FOXJ1+ Secretory" = "#CD9600", "MUC5B+" = "#8494FF", "SFTPB+ KRT5_low Basal" = "#0CB702", "SFTPB+ KRT5- Basal" = "#00A9FF", "Deuterosomal" = "#FF61CC", "SFTPC+ SCGB3A2+" = "#00B8E7") 
pdf(file.path("./", paste0("ckabnof Label Transfer to TATA colored", ".pdf")), w=8, h=8)
DimPlot(ckabnof.visual, group.by = "predicted.id", pt.size = 1, order = FALSE, cols = cols, na.value = "grey68") + NoAxes() + NoLegend()
dev.off()
pdf(file.path("./", paste0("ckdcinof Label Transfer to TATA colored", ".pdf")), w=8, h=8)
DimPlot(ckdcinof.visual, group.by = "predicted.id", pt.size = 1, order = FALSE, cols = cols, na.value = "grey68") + NoAxes() + NoLegend()
dev.off()
pdf(file.path("./", paste0("adultnof Label Transfer to TATA colored", ".pdf")), w=8, h=8)
DimPlot(adultnof.visual, group.by = "predicted.id", pt.size = 1, order = FALSE, cols = cols, na.value = "grey68") + NoAxes() + NoLegend()
dev.off()
pdf(file.path("./", paste0("Color Mapped Reference", ".pdf")), w=8, h=8)
DimPlot(epi.sub_SAE10x3.integrated_v2, group.by = "cell_type", cols = cols, pt.size = 1, order = FALSE, label = TRUE, repel = TRUE, label.size = 6) + NoAxes() + NoLegend()
dev.off()
##Data for Figure 5 G and H (graphs made in prism)
write.csv(predictions.ckdcinof.tata, "ckdcinof to tata Predictions.csv")
write.csv(predictions.ckabnof.tata, "ckabnof to tata Predictions.csv")
write.csv(predictions.adultnof.tata, "ckabnof to tata Predictions.csv")





##Merge Cells from each organoid mapping to AT2 cluster for comparison
adult.label.at2 <- subset(adultnof.mapped.tata, predicted.id == "AT2")
ckab.label.at2 <- subset(ckabnof.mapped.tata, predicted.id == "AT2")
ckdci.label.at2 <- subset(ckdcinof.mapped.tata, predicted.id == "AT2")

label.at2.merge <- merge(adult.label.at2, c(ckab.label.at2, ckdci.label.at2))
label.at2.merge.fastmnn <- label.at2.merge
set.seed(888)
DefaultAssay(label.at2.merge.fastmnn) <- "RNA"
label.at2.merge.fastmnn <- NormalizeData(label.at2.merge.fastmnn,
                                         normalization.method = "LogNormalize", scale.factor =10000)
label.at2.merge.fastmnn <-  FindVariableFeatures(label.at2.merge.fastmnn)
label.at2.merge.fastmnn <- RunFastMNN(object.list = SplitObject(label.at2.merge.fastmnn, split.by = "orig.ident"))
ElbowPlot(label.at2.merge.fastmnn, ndims = 50)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
label.at2.merge.fastmnn <- CellCycleScoring(label.at2.merge.fastmnn, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
label.at2.merge.fastmnn <- RunUMAP(label.at2.merge.fastmnn, dims = 1:30,reduction = "mnn")
label.at2.merge.fastmnn <- FindNeighbors(label.at2.merge.fastmnn, dims = 1:30,reduction = "mnn")
label.at2.merge.fastmnn <- FindClusters(label.at2.merge.fastmnn, resolution = 0.2, algorithm = 1)
label.at2.merge.fastmnn$orig.ident <- factor(label.at2.merge.fastmnn$orig.ident, levels = c("adultnof", "ckabnof", "ckdcinof"))

DimPlot(label.at2.merge.fastmnn, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 2)

##Figure 5I
pdf(file.path("./", paste0("Integrated Labeled AT2s by Sample", ".pdf")), w=11, h=8.5)
DimPlot(label.at2.merge.fastmnn, reduction = "umap",  pt.size = 2, group.by = "orig.ident")
dev.off()
#Figure 5J
pdf(file.path("./", paste0("Integrated Labeled AT2s FastMNN", ".pdf")), w=11, h=8.5)
DimPlot(label.at2.merge.fastmnn, reduction = "umap",  pt.size = 2)
dev.off()

#Figure 5K
DefaultAssay(label.at2.merge.fastmnn) <- "RNA"
label.at2.merge.fastmnn <- NormalizeData(label.at2.merge.fastmnn)
kotton.aec2.differentiation.set <- c("SFTPC", "CLDN18", "LAMP3", "SFTPB", "SFTPD", "NAPSA", "SLC34A2", "CXCL8")
kotton.aec2.maturation.set <- c("SFTPA2", "SFTPA1", "PGC", "SLPI", "CXCL5")

p <- VlnPlot(label.at2.merge.fastmnn, features = kotton.aec2.differentiation, group.by = "orig.ident",  ncol = 4, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[4]] <- p[[4]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[5]] <- p[[5]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[6]] <- p[[6]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[7]] <- p[[7]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[8]] <- p[[8]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())

pdf(file.path("./", paste0("AT2 Mapped Differentiation Kotton set", ".pdf")), w=5.7, h=1.8)
p
dev.off()

p <- VlnPlot(label.at2.merge.fastmnn, features = kotton.aec2.maturation, group.by = "orig.ident",  ncol = 4, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[4]] <- p[[4]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[5]] <- p[[5]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())

pdf(file.path("./", paste0("AT2 Mapped Maturation Kotton set", ".pdf")), w=5.7, h=1.8)
p
dev.off()





label.at2.merge.fastmnn.markers <- FindAllMarkers(label.at2.merge.fastmnn, min.pct = 0.25, logfc.threshold = 0.25)
label.at2.merge.fastmnn <- ScaleData(label.at2.merge.fastmnn)
label.at2.merge.fastmnn.markers %>% group_by(cluster)  %>% top_n(n = 20, wt = avg_log2FC) -> top10
pdf(file.path("./", paste0("Top10 Feature Heatmap v2", ".pdf")), w=22, h=17)
DoHeatmap(label.at2.merge.fastmnn, features = top10$gene) + NoLegend()
dev.off()

#Supplemental Table
zero.markers <- FindMarkers(label.at2.merge.fastmnn, ident.1 = 0, min.pct = 0.25)
one.markers <- FindMarkers(label.at2.merge.fastmnn, ident.1 = 1, min.pct = 0.25)
two.markers <- FindMarkers(label.at2.merge.fastmnn, ident.1 = 2, min.pct = 0.25)

write.csv(zero.markers, "Cluster0.csv")
write.csv(one.markers, "Cluster1.csv")
write.csv(two.markers, "Cluster2.csv")

Idents(label.at2.merge.fastmnn) <- "orig.ident"
label.at2.merge.fastmnn.markers <- FindAllMarkers(label.at2.merge.fastmnn, min.pct = 0.25, logfc.threshold = 0.25)

#Supplemental Table

write.csv(label.at2.merge.fastmnn.markers, "AT2 label transfer enrichment markers.csv")


adultnof.markers <- FindMarkers(label.at2.merge.fastmnn, ident.1 = "adultnof", min.pct = 0.25)
ckabnof.markers <- FindMarkers(label.at2.merge.fastmnn, ident.1 = "ckabnof", min.pct = 0.25)
ckdcinof.markers <- FindMarkers(label.at2.merge.fastmnn, ident.1 = "ckdcinof", min.pct = 0.25)

write.csv(adultnof.markers, "AdultEnriched.csv")
write.csv(ckabnof.markers, "ckabEnriched.csv")
write.csv(ckdcinof.markers, "ckdciEnriched.csv")

#Figure 5L
krasnow.at2.markers.2.list <- list(c("SFTPC",
                                     "NPC2",
                                     "SFTPA1",
                                     "NAPSA",
                                     "SFTPA2",
                                     "CTSH",
                                     "PGC",
                                     "SFTPD",
                                     "LAMP3",
                                     "ABCA3",
                                     "CHI3L2",
                                     "CA2",
                                     "DBI",
                                     "SERPINA1",
                                     "WIF1",
                                     "LRRK2",
                                     "C11orf96",
                                     "NRGN",
                                     "SFTA2",
                                     "PLA2G1B",
                                     "HHIP",
                                     "PEBP4",
                                     "CPB2",
                                     "NNMT",
                                     "MFSD2A",
                                     "SFTPB",
                                     "HLA-DPB1",
                                     "CRTAC1",
                                     "C2",
                                     "MALL",
                                     "MID1IP1",
                                     "SLC22A31",
                                     "FTL",
                                     "P3H2",
                                     "HLA-DPA1",
                                     "AK1",
                                     "LHFPL3-AS2",
                                     "C4BPA",
                                     "C16orf89",
                                     "CACNA2D2",
                                     "HLA-DRA",
                                     "TMEM243",
                                     "DRAM1",
                                     "LPCAT1",
                                     "TMEM163",
                                     "CXCL2",
                                     "SFTA3",
                                     "HHIP-AS1",
                                     "CD74",
                                     "ETV5",
                                     "SLC34A2",
                                     "DMBT1",
                                     "FGGY",
                                     "HLA-DRB1",
                                     "RASGRF1",
                                     "NECAB1",
                                     "SELENBP1",
                                     "SDR16C5",
                                     "TFPI",
                                     "HOPX",
                                     "RND1",
                                     "CD36",
                                     "FABP5",
                                     "MUC1",
                                     "RGS16",
                                     "ALPL",
                                     "ALOX15B",
                                     "LRRC36",
                                     "KCNJ15",
                                     "CSF3R",
                                     "SCD",
                                     "LGALSL",
                                     "LANCL1-AS1",
                                     "PPP1R1B",
                                     "SLC46A2",
                                     "DCXR",
                                     "C3",
                                     "NFKBIA",
                                     "BMP2",
                                     "DUSP6",
                                     "SLC6A14",
                                     "HLA-DRB5",
                                     "AREG",
                                     "GKN2",
                                     "CAT",
                                     "EDNRB",
                                     "CEBPD",
                                     "KCNJ8",
                                     "MSMO1",
                                     "KIAA1324L",
                                     "SNX30",
                                     "ETV1",
                                     "PARM1",
                                     "ZNF385B",
                                     "FASN",
                                     "FBP1",
                                     "HMOX1",
                                     "CITED2",
                                     "PLD3",
                                     "PMM1",
                                     "CDC42EP1",
                                     "ODC1",
                                     "ORM1",
                                     "HLA-DMA",
                                     "SPRY4",
                                     "SMAGP",
                                     "ACADL",
                                     "B3GNT8",
                                     "AGPAT2",
                                     "ESAM",
                                     "ASRGL1",
                                     "EPHX1",
                                     "LPL",
                                     "QDPR",
                                     "CISH",
                                     "MTRR",
                                     "CHI3L1",
                                     "LGMN",
                                     "CD44",
                                     "HLA-DQB1",
                                     "S100A14",
                                     "MSN",
                                     "MLPH",
                                     "GADD45B",
                                     "MBIP",
                                     "SOCS2",
                                     "GSTA4",
                                     "EP300-AS1",
                                     "TTN",
                                     "ACSL4",
                                     "ZDHHC3",
                                     "HP",
                                     "PID1",
                                     "AQP1",
                                     "HSD17B4",
                                     "SNX25",
                                     "PEBP1",
                                     "FGG",
                                     "C1orf21",
                                     "SPTSSA",
                                     "IDI1",
                                     "CXCL17",
                                     "AKAP13",
                                     "BTG1",
                                     "FMO5",
                                     "FDPS",
                                     "RAB27A",
                                     "TMSB4X",
                                     "ASAH1",
                                     "BLVRB",
                                     "FLRT3",
                                     "SECISBP2L",
                                     "CDK2AP2",
                                     "RBPMS-AS1",
                                     "TSC22D1",
                                     "MRPL14",
                                     "CHCHD7",
                                     "STC1",
                                     "ADI1",
                                     "ATP6V0E1",
                                     "NTN4",
                                     "LDHA",
                                     "TIFA",
                                     "GEM",
                                     "SAT2",
                                     "STEAP4",
                                     "SLC25A5",
                                     "TMEM41A",
                                     "POLR2C",
                                     "IFITM2",
                                     "LTA4H",
                                     "TPD52L1",
                                     "ZFP36",
                                     "ENO1",
                                     "SCP2",
                                     "CKS2",
                                     "HLA-DMB",
                                     "CYP51A1",
                                     "CHP1",
                                     "CREB3L1",
                                     "GSPT1",
                                     "CD83",
                                     "BRI3",
                                     "HMGCS1",
                                     "PTP4A3",
                                     "SOCS3",
                                     "SERPINB1",
                                     "AZGP1",
                                     "PNRC1",
                                     "SOD2",
                                     "XBP1",
                                     "GADD45G",
                                     "CSF3",
                                     "NFKBIZ",
                                     "TXNIP",
                                     "CXCL3",
                                     "PLIN2",
                                     "SEC61G",
                                     "MED24"))


epi.sub_SAE10x3.integrated_v2 <- readRDS("~/scRNAseq/Data/Cell/TATA Lab TRB-SC/epi.sub_SAE10x3.integrated_v2.RDS")
epi.sub_SAE10x3.integrated_v2 <- UpdateSeuratObject(epi.sub_SAE10x3.integrated_v2)
epi.sub_SAE10x3.integrated_v2 <- RunUMAP(object = epi.sub_SAE10x3.integrated_v2, dims = 1:20,   min.dist = 1, n.neighbors = 20, spread = 1.5, local.connectivity = 10, return.model = TRUE)
DefaultAssay(epi.sub_SAE10x3.integrated_v2) <- "RNA"




tata.at2s <- subset(epi.sub_SAE10x3.integrated_v2, idents = "AT2")
tata.mcs <- subset(epi.sub_SAE10x3.integrated_v2, idents = "Ciliated")

tata.mcs <- AddMetaData(tata.mcs, "ciliated", col.name = "orig.ident")
label.at2.merge.tata.threef.cil <- merge(tata.at2s, c(adult.label.at2, ckab.label.at2, ckdci.label.at2, threef, tata.mcs))

label.at2.merge.tata.threef.cil$orig.ident <- factor(label.at2.merge.tata.threef.cil$orig.ident, levels = c("SeuratProject","adultnof", "ckabnof", "ckdcinof", "threef", "ciliated"))


DefaultAssay(label.at2.merge.tata.threef.cil) <- "RNA"
label.at2.merge.tata.threef.cil <- NormalizeData(label.at2.merge.tata.threef.cil)

y <- label.at2.merge.tata.threef.cil
y <- AddModuleScore(y, features = krasnow.at2.markers.2.list, name = "At2KrasnowMarkers", search = TRUE, assay = "RNA")
p <- VlnPlot(y, features = "At2KrasnowMarkers1", pt.size = 0, group.by = "orig.ident")
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())

pdf(file.path("./", paste0("Alveolar Differentation Scoring 199 Krasnow Published TATA Threef Cil Included")), w=8, h=2.5)
p
dev.off()   

##Data for 5M
adultnofvckab.markers.kat2.features <- FindMarkers(label.at2.merge.fastmnn, ident.1 = "adultnof", ident.2 = "ckabnof", features = krasnow.at2.markers.2, min.pct = 0.25)
adultnofvckdci.markers.kat2.features <- FindMarkers(label.at2.merge.fastmnn, ident.1 = "adultnof", ident.2 = "ckdcinof", features = krasnow.at2.markers.2, min.pct = 0.25)

#Supplementary Figure 5g and h
DefaultAssay(epi.sub_SAE10x3.integrated_v2) <- "RNA"
epi.sub_SAE10x3.integrated_v2 <- NormalizeData(epi.sub_SAE10x3.integrated_v2)
tata.acb.goblet.markers <- FindMarkers(epi.sub_SAE10x3.integrated_v2, ident.1 = "MUC5AC+ MUC5B+", min.pct = 0.25)
tata.b.goblet.markers <- FindMarkers(epi.sub_SAE10x3.integrated_v2, ident.1 = "MUC5B+", min.pct = 0.25)
write.csv(tata.acb.goblet.markers, "TATA ACB Goblet Markers.csv")
write.csv(tata.b.goblet.markers, "TATA B Goblet Markers.csv")

ckdci.goblet <- subset(ckdcinof, idents = 2)
ckdci.goblet <- AddMetaData(ckdci.goblet, col.name = "origin", "ckdci")
ckab.goblet <- subset(ckabnof, idents = 7)
ckab.goblet <- AddMetaData(ckab.goblet, col.name = "origin", "ckab")

tata.ACB.goblet <- subset(epi.sub_SAE10x3.integrated_v2, idents = "MUC5AC+ MUC5B+")
tata.ACB.goblet <- AddMetaData(tata.ACB.goblet, col.name = "origin", "TATA ACB Goblet")
tata.B.goblet <- subset(epi.sub_SAE10x3.integrated_v2, idents = "MUC5B+")
tata.B.goblet <- AddMetaData(tata.B.goblet, col.name = "origin", "TATA B Goblet")

tata.acb.goblet.list <- list(c("BPIFB1",
                               "TFF3",
                               "LCN2",
                               "MUC5AC",
                               "SCGB3A1",
                               "MUC5B",
                               "MSMB",
                               "LTF",
                               "SCGB1A1",
                               "SLPI",
                               "SAA1",
                               "WFDC2",
                               "PRSS23",
                               "SERPINB3",
                               "CYP2F1",
                               "PIGR",
                               "TSPAN8",
                               "S100P",
                               "AGR2",
                               "CTSC",
                               "SAA2",
                               "CXCL6",
                               "TGM2",
                               "C3",
                               "TMEM45A",
                               "CXCL17",
                               "XBP1",
                               "MDK",
                               "CXCL1",
                               "VMO1",
                               "CP",
                               "CLCA2",
                               "AKR1C1",
                               "CHP2",
                               "KRT19",
                               "AQP5",
                               "SERPINA3",
                               "RARRES1",
                               "AZGP1",
                               "FAM3D",
                               "CRISP3",
                               "CXCL8",
                               "GALNT6",
                               "AKR1C2",
                               "NTS",
                               "PI3",
                               "PLEKHS1",
                               "FABP5",
                               "ADAM28",
                               "SLC4A4",
                               "CDC42EP5",
                               "PSCA",
                               "RDH10",
                               "ATP12A",
                               "RHOV",
                               "CCNO",
                               "RIMS1",
                               "KLK11",
                               "SLC9A3R2",
                               "FCGBP",
                               "IGFBP3",
                               "CD82",
                               "HES4",
                               "CPD",
                               "CFB",
                               "MUCL1",
                               "LYZ",
                               "CYBA",
                               "TNFSF10",
                               "GALNT7",
                               "PLK2",
                               "SLC31A1",
                               "ERN2",
                               "SPDEF",
                               "BACE2",
                               "ST6GAL1",
                               "FOXN3",
                               "WNK2",
                               "SLC9A3R1",
                               "KLK10",
                               "ZG16B",
                               "TENT5C",
                               "NPDC1",
                               "DUSP4",
                               "CCND1",
                               "ALDH1A3",
                               "EPAS1",
                               "MUC16",
                               "LY6D",
                               "LY6E",
                               "INSR",
                               "KRT7",
                               "SSR4",
                               "EPHX1",
                               "SLC15A2",
                               "KLF4",
                               "SULT2B1",
                               "RND3",
                               "IL6",
                               "TOX3",
                               "MMP7",
                               "LRRC26",
                               "ANPEP",
                               "WNT5A",
                               "TCEA3",
                               "NDRG2",
                               "DHRS9",
                               "CREB3L1",
                               "BMPR1B",
                               "SMIM31",
                               "CHST9",
                               "CRACR2B",
                               "MUC20",
                               "AC007681.1",
                               "KRT4",
                               "CTSB",
                               "C9orf16",
                               "GSTA1",
                               "RRBP1",
                               "GOLM1",
                               "SLC6A14",
                               "CLDN10",
                               "NRARP",
                               "LMCD1",
                               "ECE1",
                               "PAM",
                               "CAMK1D",
                               "PTGFR",
                               "PPP1R1A",
                               "FAM83D",
                               "RASSF9",
                               "SLC12A2",
                               "ADRA2A",
                               "EFNB2",
                               "TMEM205",
                               "HLA-B",
                               "CAPN13",
                               "STARD10",
                               "CRYM",
                               "SERPINB11",
                               "SORD",
                               "CEACAM6",
                               "GPX3",
                               "SPINT1",
                               "PDE8B",
                               "LIF",
                               "KLF5",
                               "TPD52L1",
                               "GCHFR",
                               "GLRX",
                               "CLDN22",
                               "CYP2J2",
                               "HLA-C",
                               "EHF",
                               "HLA-DRA",
                               "ALDH1A1",
                               "PTN",
                               "HS3ST1",
                               "PAQR4",
                               "LYN",
                               "TFPI2",
                               "HEY1",
                               "STXBP6",
                               "LINC02185",
                               "CST3",
                               "B3GALT5",
                               "GRAMD2B",
                               "TIMP1",
                               "CHL1",
                               "HCAR2",
                               "GSN",
                               "ASS1",
                               "NUDT8",
                               "ANKRD36C",
                               "ATP8B1",
                               "KLF3",
                               "CFTR",
                               "PROM1",
                               "PTPRZ1",
                               "VSIG2",
                               "EPS8L1",
                               "MAFB",
                               "SOX21",
                               "DTX4",
                               "PDLIM5",
                               "S100A9",
                               "PODXL",
                               "ST6GALNAC1",
                               "CX3CL1",
                               "MESP1",
                               "STT3B",
                               "EFHD1",
                               "KIAA1324",
                               "SERPINF1",
                               "SCNN1B",
                               "MARCKSL1",
                               "CD74",
                               "BEND5",
                               "CXCL3",
                               "FKBP11"))

tata.b.goblet.list <- list(c("SCGB3A1",
                             "BPIFB1",
                             "LCN2",
                             "LTF",
                             "SAA1",
                             "SCGB1A1",
                             "SLPI",
                             "WFDC2",
                             "C3",
                             "CXCL1",
                             "MUC5B",
                             "RARRES1",
                             "PRSS23",
                             "SERPINB3",
                             "AZGP1",
                             "SERPINA3",
                             "CYP2F1",
                             "VMO1",
                             "TGM2",
                             "TMEM45A",
                             "SAA2",
                             "CXCL8",
                             "FABP5",
                             "CRISP3",
                             "CXCL6",
                             "AKR1C1",
                             "IL6",
                             "AKR1C2",
                             "XBP1",
                             "AGR2",
                             "CP",
                             "S100A9",
                             "PIGR",
                             "MDK",
                             "TSPAN8",
                             "MUCL1",
                             "IGFBP3",
                             "KRT19",
                             "CTSC",
                             "HES4",
                             "PLEKHS1",
                             "LIF",
                             "ADAM28",
                             "GPX3",
                             "KRT7",
                             "CHP2",
                             "MMP7",
                             "CD82",
                             "CXCL17",
                             "PI3",
                             "PLAUR",
                             "ALDH1A3",
                             "RHOV",
                             "ZG16B",
                             "S100P",
                             "CFB",
                             "NTS",
                             "CXCL2",
                             "TFPI2",
                             "FCGBP",
                             "LY6D",
                             "ANPEP",
                             "CCND1",
                             "EPHX1",
                             "WNK2",
                             "HLA-DRA",
                             "IL1R1",
                             "SLC4A4",
                             "ATP12A",
                             "KRT4",
                             "TNFAIP2",
                             "SLC6A14",
                             "KLK10",
                             "PLK2",
                             "TNFSF10",
                             "NRARP",
                             "AQP5",
                             "EMP1",
                             "DUSP5",
                             "INSR",
                             "RDH10",
                             "CPD",
                             "CX3CL1",
                             "RCAN1",
                             "FOLR1",
                             "GLRX",
                             "WNT5A",
                             "CXCL3",
                             "PSCA",
                             "TOX3",
                             "ST6GAL1",
                             "EPAS1",
                             "CLCA2",
                             "FAM3D",
                             "PDE8B",
                             "BACE2",
                             "SLC9A3R2",
                             "HLA-DRB5",
                             "CD74",
                             "HLA-DRB1",
                             "ALDH1A1",
                             "SLC9A3R1",
                             "TNFAIP3",
                             "FOXN3",
                             "LYN",
                             "LY6E",
                             "EFHD1",
                             "RND3",
                             "GALNT6",
                             "CDC42EP5",
                             "SAA4",
                             "AC007681.1",
                             "GADD45A",
                             "CRABP2",
                             "TMEM205",
                             "ASS1",
                             "HS3ST1",
                             "KLF4",
                             "NDRG2",
                             "GCHFR",
                             "PAM",
                             "MSMB",
                             "CTSB",
                             "EFNB2",
                             "PTGES",
                             "MARCKSL1",
                             "KLK11",
                             "PROM1",
                             "SLC15A2",
                             "TPD52L1",
                             "SERPINB11",
                             "IRS2",
                             "ZC3H12A",
                             "PODXL",
                             "PDP1",
                             "SLC31A1",
                             "RAB5IF",
                             "SLC12A2",
                             "PTN",
                             "CAMK1D",
                             "CCNO",
                             "NPDC1",
                             "AHNAK",
                             "EDN1",
                             "SULT2B1",
                             "CYP2J2",
                             "HLA-DPA1",
                             "IQGAP2",
                             "BMPR1B",
                             "GALNT7",
                             "ERN2",
                             "SPDEF",
                             "PDLIM5",
                             "LITAF",
                             "SERPINF1",
                             "LRRC26",
                             "GRAMD2B",
                             "RHOBTB3",
                             "CFTR",
                             "CYP2A13",
                             "CNN3",
                             "RIMS1",
                             "PHLDA2",
                             "FUT2",
                             "TCEA3",
                             "SPINT1",
                             "SMIM31",
                             "HLA-C",
                             "S100A16",
                             "TENT5C",
                             "MUC16",
                             "PLAU",
                             "TRIM16",
                             "FAM3C",
                             "C9orf16",
                             "TNFRSF21",
                             "MUC20",
                             "STT3B",
                             "STEAP1",
                             "B4GALT5",
                             "PAQR4",
                             "PPP1R1A",
                             "SH3BP4",
                             "CLCF1",
                             "HLA-DPB1",
                             "ITPKC",
                             "PF4V1",
                             "NUDT8",
                             "SEC14L1",
                             "B3GALT5",
                             "CLDN22")) #194

goblet.merge <- merge(tata.ACB.goblet, c(tata.B.goblet, ckab.goblet, ckdci.goblet))
goblet.merge$origin <- factor(goblet.merge$origin, levels = c("TATA ACB Goblet", "TATA B Goblet", "ckab", "ckdci"))
y <- goblet.merge

y <- AddModuleScore(y, features = tata.acb.goblet.list, name = "acbgoblet", search = TRUE, assay = "RNA")
p <- VlnPlot(y, features = "acbgoblet1", pt.size = 0, group.by = "origin")
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank()) + NoLegend()
pdf(file.path("./", paste0("TATA ACB Goblet Score.pdf")), w=4, h=4)
p
dev.off()
y <- AddModuleScore(y, features = tata.b.goblet.list, name = "bgoblet", search = TRUE, assay = "RNA")
p <- VlnPlot(y, features = "bgoblet1", pt.size = 0, group.by = "origin")
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank()) + NoLegend()
pdf(file.path("./", paste0("TATA B Goblet Score.pdf")), w=4, h=4)
p
dev.off()

##Figure 6
ckdciday1.data <- Read10X(data.dir = "./") #BTP organoids treated with CKDCI for 1 day
ckdciday6.data <- Read10X(data.dir = "./") #BTP organoids treated with CKDCI for 6 days
ckdciday21.data <- Read10X(data.dir = "./")  #BTP organoids treated with CKDCI for 21 days

ckdciday1 <- CreateSeuratObject(counts = ckdciday1.data, project = "ckdciday1", min.cells = 3, min.features = 200)
ckdciday1 <- AddMetaData(ckdciday1, "ckdci Day 1", col.name = "group")
ckdciday1[["percent.mt"]] <- PercentageFeatureSet(ckdciday1, pattern = "^MT-")
ckdciday6 <- CreateSeuratObject(counts = ckdciday6.data, project = "ckdciday6", min.cells = 3, min.features = 200)
ckdciday6 <- AddMetaData(ckdciday6, "ckdci Day 6", col.name = "group")
ckdciday6[["percent.mt"]] <- PercentageFeatureSet(ckdciday6, pattern = "^MT-")
ckdciday21 <- CreateSeuratObject(counts = ckdciday21.data, project = "ckdciday21", min.cells = 3, min.features = 200)
ckdciday21 <- AddMetaData(ckdciday21, "ckdci Day 21", col.name = "group")
ckdciday21[["percent.mt"]] <- PercentageFeatureSet(ckdciday21, pattern = "^MT-")
ckdci <- merge(threef, c(ckdciday1, ckdciday6, ckdciday21))
ckdci.fastmnn <- subset(ckdci, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
set.seed(888)
DefaultAssay(ckdci.fastmnn) <- "RNA"
ckdci.fastmnn <- NormalizeData(ckdci.fastmnn,
                               normalization.method = "LogNormalize", scale.factor =10000)
ckdci.fastmnn <-  FindVariableFeatures(ckdci.fastmnn)
ckdci.fastmnn <- RunFastMNN(object.list = SplitObject(ckdci.fastmnn, split.by = "orig.ident"))
ElbowPlot(ckdci.fastmnn, ndims = 50)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ckdci.fastmnn <- CellCycleScoring(ckdci.fastmnn, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ckdci.fastmnn <- RunUMAP(ckdci.fastmnn, dims = 1:30,reduction = "mnn")
ckdci.fastmnn <- FindNeighbors(ckdci.fastmnn, dims = 1:30,reduction = "mnn")
ckdci.fastmnn <- FindClusters(ckdci.fastmnn, resolution = 0.3, algorithm = 1)
DimPlot(ckdci.fastmnn, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 2)



##Data for Figure 6a (graph made in prism)
PrctCellExpringGene(ckdci.fastmnn, genes = c("SFTPC", "SFTPA1", "SFTPB"), group.by = "orig.ident")


##Figure 6b
pdf(file.path("./", paste0("ckdci FastMNN by Sample", ".pdf")), w=11, h=8.5)
DimPlot(ckdci.fastmnn, reduction = "umap",  pt.size = 2, group.by = "group")
dev.off()
##Figure 6c
pdf(file.path("./", paste0("ckdci FastMNN", ".pdf")), w=11, h=8.5)
DimPlot(ckdci.fastmnn, reduction = "umap",  pt.size = 2, label = FALSE)
dev.off()
##


DefaultAssay(ckdci.fastmnn) <- "RNA"
ckdci.fastmnn <- NormalizeData(ckdci.fastmnn)

##Figure 6d
pdf(file.path("./", paste0("ckdci Lineage Markers from Fig 3J Dotplot", ".pdf")), w=12.5, h=5.75)
Dotplot_Zhiwei_Version(ckdci.fastmnn, ckab.fastmnn.markers)
dev.off()
##Supplementary Figure 6a (graph made in PRISM)
write.csv(table(ckdci.fastmnn$group, ckdci.fastmnn$Phase), "RunbyPhase.csv")
 
##Supplementary Figure 6b
ckab.fastmnn.cluster.6 <- subset(ckab.fastmnn, ident = 6)
ckdci.fastmnn.cluster.5 <- subset(ckdci.fastmnn, ident = 5)
diff.neuroendocrine.compare <- merge(ckab.fastmnn.cluster.6, ckdci.fastmnn.cluster.5)

DefaultAssay(diff.neuroendocrine.compare) <- "RNA"
diff.neuroendocrine.compare <- NormalizeData(diff.neuroendocrine.compare)

p <- VlnPlot(diff.neuroendocrine.compare, features = c("SFTPC", "ASCL1", "CHGA"),  ncol = 3, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())


pdf(file.path("./", paste0("Neuroendocrine Comparison", ".pdf")), w=5.7, h=1.8)
p
dev.off()
#statistics for Supplementary Figure 6b
diff.neroendocrine.compare.markers <- FindAllMarkers(diff.neuroendocrine.compare, min.pct = 0.25)
write.csv(diff.neroendocrine.compare.markers, "NEClusterDiffMarkers.csv")
##Supplementary Figure 6c
pdf(file.path("./", paste0("ckdci FastMNN Dotplot Cluster 3 and 6", ".pdf")), w=12.5, h=5.75)
Dotplot_Zhiwei_Version(ckdci.fastmnn, c("SFTPC", "SFTPA1", "SFTPB", "NDRG1", "SLC5A3", "SLC2A3", "SLC6A6", "CXCL8", "TF", "LCN2"))
dev.off()

pdf(file.path("./", paste0("ckdci FastMNN All Highlight", ".pdf")), w=11, h=8.5)
DimPlot(ckdci.fastmnn, cells.highlight = list(ckdciday21.bestat2.cellids, ckdciday6.bestat2.cellids, ckdciday1.bestat2.cellids, threef.cellids), cols.highlight = c("purple", "blue", "green", "yellow"), pt.size = 2, sizes.highlight = 2, order = TRUE)
dev.off()

##Figure 6e
ckdci.fastmnn <- FindClusters(ckdci.fastmnn, resolution = 1, algorithm = 2)
DimPlot(ckdci.fastmnn, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 2)

DefaultAssay(ckdci.fastmnn) <- "RNA"
ckdci.fastmnn <- NormalizeData(ckdci.fastmnn)

sce.ckdci.nosct <- as.SingleCellExperiment(ckdci.fastmnn)
sce.ckdci.nosct <- slingshot(sce.ckdci.nosct, clusterLabels = sce.ckdci.nosct@colData@listData[["ident"]]
                             , reducedDim = "UMAP",
                             start.clus = "6", allow.breaks = FALSE, stretch = 0, omega = TRUE, extend = "n")
lnes.ckdci <- getLineages(reducedDim(sce.ckdci.nosct,"UMAP"), sce.ckdci.nosct@colData@listData[["ident"]], start.clus = "6")
pt.ckdci.nosct <- slingPseudotime(sce.ckdci.nosct, na = TRUE)
pt.ckdci.nosct[is.na(pt.ckdci.nosct)] <- 0
colors <- pal[cut(pt.ckdci.nosct[,1], breaks = 100)]

pdf(file.path("./", paste0("ckdci Psuedotime with Lineage", ".pdf")), w=11, h=8.5)
plot(reducedDim(sce.ckdci.nosct, "UMAP"), col = colors, pch = 16, cex = 1)
plot(SlingshotDataSet(sce.ckdci.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 2, add = TRUE)
plot(SlingshotDataSet(sce.ckdci.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 3, add = TRUE)
plot(SlingshotDataSet(sce.ckdci.nosct), lwd =6, type = 'lineage', col = c("gray"), linInd = 5, add = TRUE)
plot(SlingshotDataSet(sce.ckdci.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 4, add = TRUE)
plot(SlingshotDataSet(sce.ckdci.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 6, add = TRUE)
plot(SlingshotDataSet(sce.ckdci.nosct), lwd=8, type = 'lineage', col = c("red"), linInd = 1, add = TRUE)
dev.off()

##Figure 6f
pdf(file.path("./", paste0("ckdci FastMNN SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdci.fastmnn, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()


pdf(file.path("./", paste0("ckdci FastMNN SCGB3A2", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdci.fastmnn, features = "SCGB3A2", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

pdf(file.path("./", paste0("ckdci FastMNN RNASE1", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdci.fastmnn, features = "RNASE1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

##Figure 6g
alvdiff.fastmnn$group <- factor(alvdiff.fastmnn$group, c("2FAB Day 21", "2FAB Day 6", "2FAB Day 1", "3F", "ckdci Day 1", "ckdci Day 6", "ckdci Day 21"))
p <- VlnPlot(alvdiff.fastmnn, features = "SCGB3A2", group.by = "group",  pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN SCGB3A2 set", ".pdf")), w=8, h=4)
p
dev.off()
p <- VlnPlot(alvdiff.fastmnn, features = "SFTPC", group.by = "group",  pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN SFTPC set", ".pdf")), w=6, h=4)
p
dev.off()
p <- VlnPlot(alvdiff.fastmnn, features = "RNASE1", group.by = "group",  pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN RNASE1 set", ".pdf")), w=6, h=4)
p
dev.off()

##Figure 6h
alvdiff.combined <- merge(threef, c(ckabday1, ckabday6, ckabday21, ckdciday1, ckdciday6, ckdciday21))
alvdiff.combined <- subset(alvdiff.combined, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
DefaultAssay(alvdiff.combined) <- "RNA"
alvdiff.combined <- NormalizeData(alvdiff.combined)
p <- VlnPlot(alvdiff.combined, features = "LAMP3", group.by = "group",  pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN LAMP3 set", ".pdf")), w=6, h=4)
p
dev.off()
p <- VlnPlot(alvdiff.combined, features = "SFTPB", group.by = "group",  pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN LAMP3 set", ".pdf")), w=6, h=4)
p
dev.off()
p <- VlnPlot(alvdiff.combined, features = "SFTPA1", group.by = "group",  pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN SFTPA1 set", ".pdf")), w=6, h=4)
p
dev.off()


##Figure 6j - m
###Label transfer cells from all timepoints of CK + AB and CK+DCI differentiations to same reference dataset as used in figure 5 (Murthy et al., 2022)
epi.sub_SAE10x3.integrated_v2 <- readRDS("~/scRNAseq/Data/Cell/TATA Lab TRB-SC/epi.sub_SAE10x3.integrated_v2.RDS")
epi.sub_SAE10x3.integrated_v2 <- UpdateSeuratObject(epi.sub_SAE10x3.integrated_v2)
epi.sub_SAE10x3.integrated_v2 <- RunUMAP(object = epi.sub_SAE10x3.integrated_v2, dims = 1:20,   min.dist = 1, n.neighbors = 20, spread = 1.5, local.connectivity = 10, return.model = TRUE)
DefaultAssay(epi.sub_SAE10x3.integrated_v2) <- "RNA"


ckdci <- merge(ckdciday1, c(ckdciday6, ckdciday21))

ckdci <- subset(ckdci, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
ckdci <- UpdateSeuratObject(ckdci)

ckdci <- NormalizeData(ckdci)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ckdci <- CellCycleScoring(ckdci, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ckdci <- FindVariableFeatures(ckdci, selection.method = "vst", nfeatures = 2000)
ckdci <- ScaleData(ckdci)
ckdci <- RunPCA(ckdci, features = VariableFeatures(ckdci))
ElbowPlot(ckdci, ndims = 50)
ckdci <- RunUMAP(ckdci, dims = 1:20, reduction.name = "umap", return.model = TRUE)
ckdci <- FindNeighbors(ckdci, dims = 1:20)
ckdci <- FindClusters(ckdci, resolution = 0.3)
DimPlot(ckdci, label = TRUE, pt.size = 1)

common.features <- intersect(rownames(epi.sub_SAE10x3.integrated_v2), rownames(ckdci))
common.features.2 <- intersect(rownames(common.features, rownames(ckab)))
length(x = common.features.2)
epi.sub_SAE10x3.integrated_v2.anchors <- FindTransferAnchors(
  reference = epi.sub_SAE10x3.integrated_v2,
  normalization.method = "LogNormalize",
  query = ckdci,
  reference.reduction = "pca",
  features = common.features.2,
  dims = 1:20
)


ckdci.mapped.tata <- MapQuery(
  anchorset = epi.sub_SAE10x3.integrated_v2.anchors,
  query = ckdci,
  reference = epi.sub_SAE10x3.integrated_v2,
  refdata = "cell_type",
  reference.reduction = "pca", 
  transferdata.args = list(k.weight = 10),
  integrateembeddings.args = list(k.weight = 10),
  reduction.model = "umap"
)



ckdci.mapped.tata@reductions[["umap"]] <- ckdci.mapped.tata@reductions[["ref.umap"]]
ckdci.visual <- merge(epi.sub_SAE10x3.integrated_v2, ckdci.mapped.tata, merge.dr = "umap")




epi.sub_SAE10x3.integrated_v2 <- readRDS("~/scRNAseq/Data/Cell/TATA Lab TRB-SC/epi.sub_SAE10x3.integrated_v2.RDS")
epi.sub_SAE10x3.integrated_v2 <- UpdateSeuratObject(epi.sub_SAE10x3.integrated_v2)
epi.sub_SAE10x3.integrated_v2 <- RunUMAP(object = epi.sub_SAE10x3.integrated_v2, dims = 1:20,   min.dist = 1, n.neighbors = 20, spread = 1.5, local.connectivity = 10, return.model = TRUE)
DefaultAssay(epi.sub_SAE10x3.integrated_v2) <- "RNA"


ckab <- merge(ckabday1, c(ckabday6, ckabday21))

ckab <- subset(ckab, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
ckab <- UpdateSeuratObject(ckab)

ckab <- NormalizeData(ckab)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ckab <- CellCycleScoring(ckab, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ckab <- FindVariableFeatures(ckab, selection.method = "vst", nfeatures = 2000)
ckab <- ScaleData(ckab)
ckab <- RunPCA(ckab, features = VariableFeatures(ckab))
ElbowPlot(ckab, ndims = 50)
ckab <- RunUMAP(ckab, dims = 1:20, reduction.name = "umap", return.model = TRUE)
ckab <- FindNeighbors(ckab, dims = 1:20)
ckab <- FindClusters(ckab, resolution = 0.3)
DimPlot(ckab, label = TRUE, pt.size = 1)


epi.sub_SAE10x3.integrated_v2.anchors <- FindTransferAnchors(
  reference = epi.sub_SAE10x3.integrated_v2,
  normalization.method = "LogNormalize",
  query = ckab,
  reference.reduction = "pca",
  features = common.features.2,
  dims = 1:20
)


ckab.mapped.tata <- MapQuery(
  anchorset = epi.sub_SAE10x3.integrated_v2.anchors,
  query = ckab,
  reference = epi.sub_SAE10x3.integrated_v2,
  refdata = "cell_type",
  reference.reduction = "pca", 
  transferdata.args = list(k.weight = 10),
  integrateembeddings.args = list(k.weight = 10),
  reduction.model = "umap"
)

ckab.mapped.tata@reductions[["umap"]] <- ckab.mapped.tata@reductions[["ref.umap"]]
ckab.visual <- merge(epi.sub_SAE10x3.integrated_v2, ckab.mapped.tata, merge.dr = "umap")

#Figure 6j
cols <- c("AT1" = "#f5fc86", "AT2" = "#7CAE00", "Proliferating" = "#ED68ED", "Proliferating AT2" ="#FF61CC", "Differentiating Basal" = "#8494FF", "Neuroendocrine" = "#ab467a" , "SFTPB+ SCGB3A2+ SCGB1A1-" = "#ABA300", "SFTPB+ SCGB3A2+ SCGB1A1+" = "#E68613", "SFTPB- KRT5+ Basal" = "#00C19A",  "Proliferating" = "#ED68ED", "MUC5AC+ MUC5B+" = "#00BE67", "Immature AT1" ="#0CB702", "Ciliated" ="#F8766D", "FOXJ1+ Secretory" = "#CD9600", "MUC5B+" = "#8494FF", "SFTPB+ KRT5_low Basal" = "#0CB702", "SFTPB+ KRT5- Basal" = "#00A9FF", "Deuterosomal" = "#FF61CC", "SFTPC+ SCGB3A2+" = "#00B8E7") 


pdf(file.path("./", paste0("ckdci Label Transfer to TATA colored", ".pdf")), w=8, h=8)
DimPlot(ckdci.visual, group.by = "predicted.id", pt.size = 1, order = FALSE, cols = cols, na.value = "grey68") + NoAxes() + NoLegend()
dev.off()
pdf(file.path("./", paste0("Color Mapped Reference", ".pdf")), w=8, h=8)
DimPlot(epi.sub_SAE10x3.integrated_v2, group.by = "cell_type", cols = cols, pt.size =1, order = FALSE, label = TRUE, label.size = 6, repel = TRUE) + NoAxes() + NoLegend()
dev.off()

pdf(file.path("./", paste0("ckab Label Transfer to TATA colored", ".pdf")), w=8, h=8)
DimPlot(ckab.visual, group.by = "predicted.id", pt.size = 1, order = FALSE, cols = cols, na.value = "grey68") + NoAxes() + NoLegend()
dev.off()




##Figure 6j


ckab.mapped.tata$group<- factor(ckab.mapped.tata$group, levels = c("Day 1", "Day 6", "Day 21")) 

ckdci.mapped.tata$orig.ident<- factor(ckdci.mapped.tata$group, levels = c("Day 1", "Day 6", "Day 21"))
pdf(file.path("./", paste0("ckab LT Assembly", ".pdf")), w=8, h=8)
DimPlot(ckab.mapped.tata, reduction = "ref.umap", group.by = "predicted.id", label = FALSE, na.value = "grey68",
        label.size = 3, repel = TRUE, pt.size = 1) + NoLegend() + NoAxes()
dev.off()


pdf(file.path("./", paste0("ckdci LT Assembly", ".pdf")), w=8, h=8)
DimPlot(ckdci.mapped.tata, reduction = "ref.umap", group.by = "predicted.id", label = FALSE, na.value = "grey68",
              label.size = 3, repel = TRUE, pt.size = 1) + NoLegend() + NoAxes()
dev.off()
##Data for Figure 6k and l
predictions.ckab.tata <- table(ckab.mapped.tata$predicted.id)
write.csv(predictions.ckab.tata, "ckab to tata Predictions.csv")
predictions.ckdci.tata <- table(ckdci.mapped.tata$predicted.id)
write.csv(predictions.ckdci.tata, "ckdci to tata Predictions.csv")

#Figure 6m
pdf(file.path("./", paste0("ckab LT Assembly", ".pdf")), w=8, h=8)
DimPlot(ckab.visual,  group.by = "group", cols = c("Day 1" = "#F8766D", "Day 6" = "#00BA38", "Day 21" = "#619CFF"), label = FALSE,
        na.value = "grey68", label.size = 3, repel = TRUE, pt.size = 1)  + NoLegend()+ NoAxes()
dev.off()
pdf(file.path("./", paste0("ckdci LT Assembly", ".pdf")), w=8, h=8)
DimPlot(ckdci.visual,  group.by = "group", cols = c("ckdciday1" = "#F8766D", "ckdciday6" = "#00BA38", "ckdciday21" = "#619CFF"), label = FALSE,
na.value = "grey68", label.size = 3, repel = TRUE, pt.size = 1)  + NoLegend()+ NoAxes()
dev.off()








##Supplementary Figure 6 d - f
##Analyze each timepoint of CKDCI differentiation separately
###Day 1
ckdciday1 <- CreateSeuratObject(counts = ckdciday1.data, project = "ckdciday1", min.cells = 3, min.features = 200)
ckdciday1 <- AddMetaData(ckdciday1, "ckdci Day 1", col.name = "group")
ckdciday1[["percent.mt"]] <- PercentageFeatureSet(ckdciday1, pattern = "^MT-")
ckdciday1 <- RenameCells(object = ckdciday1, add.cell.id = "ckdciDAY1")
VlnPlot(ckdciday1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ckdciday1 <- subset(ckdciday1, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ckdciday1 <- CellCycleScoring(ckdciday1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ckdciday1 <- SCTransform(ckdciday1, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
ckdciday1 <- RunPCA(ckdciday1, features = VariableFeatures(object = ckdciday1))
ElbowPlot(ckdciday1, ndims = 50)
DefaultAssay(ckdciday1) <- "SCT"
ckdciday1 <- RunUMAP(ckdciday1, dims = 1:18)
ckdciday1 <- FindNeighbors(ckdciday1, dims = 1:18)
ckdciday1 <- FindClusters(ckdciday1, resolution = 0.5)
DimPlot(ckdciday1, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)

##Supplementary Figure 6d
pdf(file.path("./", paste0("Day 1 Louvain", ".pdf")), w=11, h=8.5)
DimPlot(ckdciday1, pt.size = 2, label = FALSE)
dev.off()
DefaultAssay(ckdciday1) <- "RNA"
ckdciday1 <- NormalizeData(ckdciday1)

pdf(file.path("./", paste0("Day 1 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdciday1, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 1 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdciday1, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
ckdciday1.bestat2 <- subset(ckdciday1, idents = 5)

pdf(file.path("./", paste0("Day 1 Louvain Best AT2s", ".pdf")), w=11, h=8.5)
DimPlot(ckdciday1, cells.highlight = ckdciday1.bestat2.cellids, pt.size = 2, sizes.highlight = 2, cols.highlight = "lightblue")
dev.off()

###Day 6
ckdciday6 <- CreateSeuratObject(counts = ckdciday6.data, project = "ckdciday6", min.cells = 3, min.features = 200)
ckdciday6 <- AddMetaData(ckdciday6, "ckdci Day 6", col.name = "group")
ckdciday6[["percent.mt"]] <- PercentageFeatureSet(ckdciday6, pattern = "^MT-")
ckdciday6 <- RenameCells(object = ckdciday6, add.cell.id = "ckdciDAY6")
VlnPlot(ckdciday6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ckdciday6 <- subset(ckdciday6, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ckdciday6 <- CellCycleScoring(ckdciday6, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ckdciday6 <- SCTransform(ckdciday6, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
ckdciday6 <- RunPCA(ckdciday6, features = VariableFeatures(object = ckdciday6))
ElbowPlot(ckdciday6, ndims = 50)
DefaultAssay(ckdciday6) <- "SCT"
ckdciday6 <- RunUMAP(ckdciday6, dims = 1:18)
ckdciday6 <- FindNeighbors(ckdciday6, dims = 1:18)
ckdciday6 <- FindClusters(ckdciday6, resolution = 0.5)
DimPlot(ckdciday6, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)

##Supplementary Figure 6e
pdf(file.path("./", paste0("Day 6 Louvain", ".pdf")), w=11, h=8.5)
DimPlot(ckdciday6, pt.size = 2, label = FALSE)
dev.off()
DefaultAssay(ckdciday6) <- "RNA"
ckdciday6 <- NormalizeData(ckdciday6)

pdf(file.path("./", paste0("Day 6 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdciday6, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 6 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdciday6, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
ckdciday6.bestat2 <- subset(ckdciday6, idents = 1)

pdf(file.path("./", paste0("Day 6 Louvain Best AT2s", ".pdf")), w=11, h=8.5)
DimPlot(ckdciday6, cells.highlight = ckdciday6.bestat2.cellids, pt.size = 2, sizes.highlight = 2, cols.highlight = "lightblue")
dev.off()

##Day 21
ckdciday21 <- CreateSeuratObject(counts = ckdciday21.data, project = "ckdciday21", min.cells = 3, min.features = 200)
ckdciday21 <- AddMetaData(ckdciday21, "ckdci Day 21", col.name = "group")
ckdciday21[["percent.mt"]] <- PercentageFeatureSet(ckdciday21, pattern = "^MT-")
ckdciday21 <- RenameCells(object = ckdciday21, add.cell.id = "ckdciDAY21")
VlnPlot(ckdciday21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ckdciday21 <- subset(ckdciday21, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ckdciday21 <- CellCycleScoring(ckdciday21, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ckdciday21 <- SCTransform(ckdciday21, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
ckdciday21 <- RunPCA(ckdciday21, features = VariableFeatures(object = ckdciday21))
ElbowPlot(ckdciday21, ndims = 50)
DefaultAssay(ckdciday21) <- "SCT"
ckdciday21 <- RunUMAP(ckdciday21, dims = 1:18)
ckdciday21 <- FindNeighbors(ckdciday21, dims = 1:18)
ckdciday21 <- FindClusters(ckdciday21, resolution = 0.5)
DimPlot(ckdciday21, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)

##Supplementary Figure 6f

pdf(file.path("./", paste0("Day 21 Louvain", ".pdf")), w=11, h=8.5)
DimPlot(ckdciday21, pt.size = 2, label = FALSE)
dev.off()
DefaultAssay(ckdciday21) <- "RNA"
ckdciday21 <- NormalizeData(ckdciday21)

pdf(file.path("./", paste0("Day 21 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdciday21, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 21 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(ckdciday21, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
ckdciday21.bestat2 <- subset(ckdciday21, idents = 0)

pdf(file.path("./", paste0("Day 21 Louvain Best AT2s", ".pdf")), w=11, h=8.5)
DimPlot(ckdciday21, cells.highlight = ckdciday21.bestat2.cellids, pt.size = 2, sizes.highlight = 2, cols.highlight = "lightblue")
dev.off()




##Supplementary Figure 6g

ckdciday1.bestat2.cellids <- WhichCells(ckdciday1, idents = 5)
ckdciday6.bestat2.cellids <- WhichCells(ckdciday6, idents = 1)
ckdciday21.bestat2.cellids <- WhichCells(ckdciday21, idents = 0)


pdf(file.path("./", paste0("ckdci FastMNN All Highlight", ".pdf")), w=11, h=8.5)
DimPlot(ckdci.fastmnn, cells.highlight = list(ckdciday21.bestat2.cellids, ckdciday6.bestat2.cellids, ckdciday1.bestat2.cellids, threef.cellids), cols.highlight = c("purple", "blue", "green", "yellow"), pt.size = 2, sizes.highlight = 2, order = TRUE)
dev.off()
