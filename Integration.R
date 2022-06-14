library(Seurat)
library(harmony)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(MAST)
library(DoubletFinder)
library(SeuratWrappers)
library(monocle3)
library(patchwork)



########################################
########################################
files <- list.files(path = "~/Analysis/Mouse Cell Atlas/5435866/rmbatch_dge", full.name = TRUE)
anns = read.csv("MCA_CellAssignments.csv")
rownames(anns) = paste(anns[,2])

list = list()
for (i in 1: length(files)){
  
  h <- files[i]
  data <- read.table(h)
  data <- data.frame(data)

  list[i] <- CreateSeuratObject(counts = data, project = "MCA")
}
MCA = list[[1]]
for (i in 2: length(files)){
	MCA <- merge(MCA, y = list[[i]], project = "MCA")
}

saveRDS(MCA, file = "preMCA.rds")

########################################
########################################
MCA <- NormalizeData(MCA, normalization.method = "LogNormalize", scale.factor = 10000)
MCA <- FindVariableFeatures(MCA)

MCA[["percent.mt"]] <- PercentageFeatureSet(MCA, pattern = "^mt-")
MCA <- ScaleData(MCA, vars.to.regress = "percent.mt")
MCA <- RunPCA(MCA, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
MCA <- FindNeighbors(MCA, reduction = "pca", dims = 1:75, nn.eps = 0.5)
MCA <- FindClusters(MCA, resolution = 3, n.start = 10)
MCA <- RunUMAP(MCA, dims = 1:75, min.dist = 0.75)
MCA@meta.data <- anns

p <- DimPlot(MCA, reduction = "umap", pt.size = 0.1, group.by = "orig.ident", label = TRUE) + ggtitle(label = "UMAP")
p <- AugmentPlot(plot = p)
p + NoLegend()

saveRDS(MCA, file = "MCA.rds")

unique(MCA@meta.data$Tissue)
# [1] Bladder                        Bone-Marrow_c-kit              Bone_Marrow_Mesenchyme        
# [4] Bone-Marrow                    Brain                          Embryonic-Mesenchyme          
# [7] Embryonic-Stem-Cell            Female_Fetal_Gonad             Fetal_Brain                   
#[10] Fetal_Intestine                Fetal-Liver                    Fetal_Lung                    
#[13] Fetal_Stomache                 Kidney                         Liver                         
#[16] Lung                           Male_Fetal_Gonad               Mesenchymal-Stem-Cell-Cultured
#[19] Muscle                         Neonatal_Brain                 Neonatal-Calvaria             
#[22] Neonatal-Heart                 Neonatal-Muscle                Neonatal-Rib                  
#[25] Neonatal-Skin                  Ovary                          Pancreas                      
#[28] Peripheral_Blood               Placenta                       Prostate                      
#[31] Small-Intestine                Spleen                         Stomach                       
#[34] Testis                         Thymus                         Trophoblast-Stem-Cell         
#[37] Uterus                         Fetal_Kidney                   MammaryGland.Involution       
#[40] MammaryGland.Lactation         MammaryGland.Pregnancy         MammaryGland.Virgin           
#[43] NeonatalPancreas              


BREAST = subset(MCA, subset = Tissue %in% c("MammaryGland.Involution", "MammaryGland.Lactation", "MammaryGland.Pregnancy", "MammaryGland.Virgin"))
BREAST <- RunUMAP(BREAST, dims = 1:75, min.dist = 0.75)
BREAST <- FindNeighbors(BREAST, reduction = "pca", dims = 1:75, nn.eps = 0.5)
BREAST <- FindClusters(BREAST, resolution = 3, n.start = 10)

DimPlot(BREAST, group.by = "Tissue", label = TRUE) + NoLegend()
DimPlot(BREAST, group.by = "Annotation", label = TRUE) + NoLegend()

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
BREAST <- CellCycleScoring(BREAST, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Idents(object = BREAST) <- BREAST@meta.data$Annotation

DimPlot(BREAST, group.by = "Phase")

saveRDS(BREAST, file = "BREAST.rds")


modmeta <- BREAST@meta.data
IDENT <- as.character(BREAST@meta.data$Annotation)

f = read.csv("./Cellclass/EPI/Dividing Cell.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "DividingEpithelial"
f = read.csv("./Cellclass/EPI/Ductal luminal.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "DuctalLuminal"
f = read.csv("./Cellclass/EPI/Luminal progenitor.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "LuminalProgenitor"
f = read.csv("./Cellclass/EPI/Luminal.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "Luminal"
f = read.csv("./Cellclass/EPI/Myoepithelial.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "BasalMyoepithelial"
f = read.csv("./Cellclass/EPI/Secretory Alveoli.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "SecretoryAlveoli"
f = read.csv("./Cellclass/EPI/Stem and Progenitor.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "Stem&Progenitor"
f = read.csv("./Cellclass/Stromal/Adipocyte.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "Adipocyte"
f = read.csv("./Cellclass/Stromal/Bcell.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "Bcelll"
f = read.csv("./Cellclass/Stromal/Dendritic.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "Dendriticcell"
f = read.csv("./Cellclass/Stromal/Endothelilal.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "Endothelilal"
f = read.csv("./Cellclass/Stromal/Macrophage.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "Macrophage"
f = read.csv("./Cellclass/Stromal/Muscle.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "Muscle"
f = read.csv("./Cellclass/Stromal/Myeloid.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "Myeloid"
f = read.csv("./Cellclass/Stromal/NK.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "NK"
f = read.csv("./Cellclass/Stromal/Stromal.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "Stromal"
f = read.csv("./Cellclass/Stromal/Tcell.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "Tcell"
f = read.csv("./Cellclass/Stromal/Neuron.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "Neuron"
f = read.csv("./Cellclass/Stromal/Mast cell.csv")
clus = unlist(f)
IDENT[IDENT %in% clus] <- "MAST"

modmeta <- cbind(modmeta, CLASS = IDENT)
BREAST@meta.data <- modmeta


DimPlot(BREAST, group.by = "CLASS", label = TRUE) + NoLegend()

saveRDS(BREAST, file = "BREAST.rds")




VlnPlot(BREAST, "Col1a1", group.by = "CLASS", split.by = "Tissue")






COL = subset(BREAST, subset = CLASS %in% c("Bcelll", "Macrophage", "Tcell", "Stromal", "NK", "Dendriticcell", "Myeloid"))
COL <- RunUMAP(COL, dims = 1:75, min.dist = 0.75)
COL <- FindNeighbors(COL, reduction = "pca", dims = 1:75, nn.eps = 0.5)
COL <- FindClusters(COL, resolution = 3, n.start = 10)

DimPlot(COL, group.by = "CLASS", label = TRUE) + NoLegend()


saveRDS(COL, file = "COL.rds")


