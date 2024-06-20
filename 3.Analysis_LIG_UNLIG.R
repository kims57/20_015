### PART 3: Analysis

library(Seurat)
library(ggplot2)
library(clustree)
library(patchwork)
library(openxlsx)

setwd("~/20_015/data/")


## Load Datasets & Group by Tx Time & Percent Mitochondrial Filtering 
Day0 <- readRDS("~/20_015/data/cleaned_Day0.RDS")
Day0$tx_time <- "Day0"
Day0 <- PercentageFeatureSet(Day0, pattern = "^mt-", col.name = "percent.mt")
#VlnPlot(Day0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Day0 <- subset(Day0, subset = nFeature_RNA >= 200 & percent.mt < 10)

Day0_UNL <- readRDS("~/20_015/data/cleaned_Day0_UNL.RDS")
Day0_UNL$tx_time <- "Day0"
Day0_UNL <- PercentageFeatureSet(Day0_UNL, pattern = "^mt-", col.name = "percent.mt")
#VlnPlot(Day0_UNL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Day0_UNL <- subset(Day0_UNL, subset = nFeature_RNA >= 200 & percent.mt < 7)

Day0p5 <- readRDS("~/20_015/data/cleaned_Day0.5.RDS")
Day0p5$tx_time <- "Day0p5"
Day0p5 <- PercentageFeatureSet(Day0p5, pattern = "^mt-", col.name = "percent.mt")
#VlnPlot(Day0p5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Day0p5 <- subset(Day0p5, subset = nFeature_RNA >= 200 & percent.mt < 8)

Day1 <- readRDS("~/20_015/data/cleaned_Day1.RDS")
Day1$tx_time <- "Day1"
Day1 <- PercentageFeatureSet(Day1, pattern = "^mt-", col.name = "percent.mt")
#VlnPlot(Day1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Day1 <- subset(Day1, subset = nFeature_RNA >= 200 & percent.mt < 7)

Day1_LIG <- readRDS("~/20_015/data/cleaned_Day1_LIG.RDS")
Day1_LIG$tx_time <- "Day1"
Day1_LIG <- PercentageFeatureSet(Day1_LIG, pattern = "^mt-", col.name = "percent.mt")
#VlnPlot(Day1_LIG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Day1_LIG <- subset(Day1_LIG, subset = nFeature_RNA >= 200 & percent.mt < 7)

Day4 <- readRDS("~/20_015/data/cleaned_Day4.RDS")
Day4$tx_time <- "Day4"
Day4 <- PercentageFeatureSet(Day4, pattern = "^mt-", col.name = "percent.mt")
#VlnPlot(Day4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Day4 <- subset(Day4, subset = nFeature_RNA >= 200 & percent.mt < 5)

Day7 <- readRDS("~/20_015/data/cleaned_Day7.RDS")
Day7$tx_time <- "Day7"
Day7 <- PercentageFeatureSet(Day7, pattern = "^mt-", col.name = "percent.mt")
#VlnPlot(Day7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Day7 <- subset(Day7, subset = nFeature_RNA >= 200 & percent.mt < 5)


## Straight Merge [file output]
all <- merge(x= Day0 , y = c(Day0_UNL, Day0p5, Day1, Day1_LIG, Day4, Day7), add.cell.ids = c("Day0","Day0_UNL","Day0p5","Day1","Day1_LIG","Day4","Day7"), project = "all")
saveRDS(all, file = "~/20_015/outputs/straight_merged_all.rds")
DefaultAssay(all)


## SC Transformation
all <- SCTransform(all, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = TRUE)


## Linear Dimensional Reduction (PCA)
all <- RunPCA(all, features = VariableFeatures(object = all))
DimPlot(all, reduction = "pca", shuffle = TRUE)
ElbowPlot(all, ndims=50)


##Optimal clustering analysis -manually done-
all <- FindNeighbors(all, dims = 1:50)
all <- FindClusters(all, resolution = c(0,0.1,0.2,0.3,0.4,0.5))
all <- RunUMAP(all, dims = 1:50)

pdf('~/20_015/outputs/all_clustree_seurat.pdf', width = 30, height = 20)
clustree(all, node_size_range=c(10,20), node_text_size = 8)
dev.off()

Idents(all)<-"SCT_snn_res.0.1" 
#all <- FindClusters(all, resolution=0.1, cluster.name = "seurat_clusters") #alternative for line above


      

## Differentially Expressed Features [File Output]
all_DE_markers <- FindAllMarkers(all, logfc.threshold = 0.5,min.pct = 0.25)
write.csv(all_DE_markers,"~/20_015/outputs/all_DE_markers.csv", row.names = FALSE)




##Cell Annotation Here -manual or programatic-




##Assigning Cell Identity to Clusters
new.cluster.ids <- c("Fibroblast", "Endothelial", "Monocyte_MacroPhage_DC", "B_Plasma", "Epithelial_Mucin", "RBC", "Lymphatic_Endothelial_Cell", "T ", "Epithelial", "Neutrophil", "b_plasma", "Neural", "Smooth_Muscle", "epithelial_mucin", "Skeletal_Muscle")
names(new.cluster.ids) <- levels(all)
all <- RenameIdents(all, new.cluster.ids)

all$celltype <- Idents(all)



##Save Original Project Prior to Edit
saveRDS(all, file = "~/20_015/outputs/backup_all_RBC.RDS")


## Remove RBC Cluster
all <- subset(all, idents = c("RBC"), invert=TRUE)


# Cell Selector Debris Removal
DimPlot(all, reduction="umap", group.by = "celltype", shuffle = TRUE)

p <- DimPlot(all, reduction="umap", group.by = "celltype", shuffle = TRUE)
selected_cells<-CellSelector(p)

all <- subset(all, cells =  selected_cells, invert=TRUE)



##ggplot [Output Files]
celltype <- DimPlot(all, reduction="umap", group.by = "celltype", shuffle = TRUE) + xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("Cell Type")
celltype
pdf(celltype, file="~/20_015/outputs/Fig1_UMAP_Clusters.pdf")
dev.off()

txtime <-DimPlot(all, reduction = "umap", group.by = "tx_time", shuffle = TRUE)+ xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("Txt Time")
txtime
pdf(txtime, file="~/20_015/outputs/Fig2_UMAP_Txt_Time.pdf", width = 10, height = 8)
dev.off()

splittxtime <- DimPlot(all, reduction = "umap", split.by = "tx_time", shuffle = TRUE)+ xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("Split Txt Time")
splittxtime
pdf(splittxtime, file="~/20_015/outputs/Fig3_UMAP_Txt_Time_Split.pdf", width = 10, height = 8)
dev.off()


#Cell Quantification: 
table(Idents(all))

prop.table(table, 1)
prop.table(mytable, 2)

pdf('~/20_015/outputs/proportions.pdf', width = 30, height = 20)
mytable <- table(all$celltype, all$tx_time)
mydf <- as.data.frame(prop.table(mytable, 2))
ggplot(mydf, aes(x=Var2, y = Freq, fill=Var1)) + 
  geom_bar(stat="identity") +
  scale_fill_discrete(name = "Cluster") +
  xlab("Time") + 
  ylab("Proportion")
dev.off()

#mytable <- table(all$celltype, all$orig.ident)
#mydf <- as.data.frame(prop.table(mytable, 2))
#ggplot(mydf, aes(x=Var2, y = Freq, fill=Var1)) + 
  #geom_bar(stat="identity") +
  #scale_fill_discrete(name = "Cluster") +
  #xlab("Time") + 
  #ylab("Proportion")




#DE Feature DotPlot: -use all_DE_markers.csv to select features-
DotPlot(all, features = c("Col1a1", "Col1a2", "Col5a2", "Cd68", "Lyz2", "Fcer1g", "S100a8", "S100a9", "Il1b", "Igkc", "H2-Aa", "Ighm", "Trbc2", "Cd3g"), group.by = "celltype") + RotatedAxis()

DotPlot(all, features = c(""), group.by = "celltype") + RotatedAxis()



## Create artificial ident column
all$celltype.timepoint <- paste(all$celltype,all$tx_time, sep = "_")
Idents(all)<-"celltype.timepoint"


sample <- c("Fibroblast", "Endothelial", "Monocyte_MacroPhage_DC", "B_Plasma", "Epithelial_Mucin", "RBC", "Lymphatic_Endothelial_Cell", "T ", "Epithelial", "Neutrophil", "b_plasma", "Neural", "Smooth_Muscle", "epithelial_mucin", "Skeletal_Muscle")
for(x in sample) {
  p5_DE_genes <- FindMarkers(all, ident.1 = c(paste0(x, "_Day0")), ident.2 = c(paste0(x, "_Day0p5")), test.use="MAST")
  print("THE DE portion finished succesfully")
  write.csv(p5_DE_genes, paste0(x,"_p5_early_DE_markers.csv"), row.names = FALSE)
  
  one_DE_genes <- FindMarkers(all, ident.1 = c(paste0(x, "_Day0")), ident.2 = c( paste0(x, "_Day1")), test.use="MAST")
  write.csv(one_DE_genes, paste0(x,"_1_DE_markers.csv"), row.names = FALSE)
  
  four_DE_genes <- FindMarkers(all, ident.1 = c(paste0(x, "_Day0")), ident.2 = c( paste0(x, "_Day4")), test.use="MAST")
  write.csv(four_DE_genes, paste0(x,"_4_DE_markers.csv"), row.names = FALSE)
  
  seven_DE_genes <- FindMarkers(all, ident.1 = c(paste0(x, "_Day0")), ident.2 =  c(paste0(x, "_Day7")), test.use="MAST")
  write.csv(seven_DE_genes, paste0(x,"_7_DE_markers.csv"), row.names = FALSE)
  
  rm(p5_DE_genes, one_DE_genes, four_DE_genes, seven_DE_genes)
  gc()
}


Idents(all)

#Cluster 0
Fibroblast_p5_DE_genes <- FindMarkers(all, ident.1 = c("Fibroblast_Day0"), ident.2 = c("Fibroblast_Day0p5"), test.use="MAST")
Fibroblast_1_DE_genes <- FindMarkers(all, ident.1 = c("Fibroblast_Day0"), ident.2 = c("Fibroblast_Day1"), test.use="MAST")
Fibroblast_4_DE_genes <- FindMarkers(all, ident.1 = c("Fibroblast_Day0"), ident.2 = c("Fibroblast_Day4"), test.use="MAST")
Fibroblast_7_DE_genes <- FindMarkers(all, ident.1 = c("Fibroblast_Day0"), ident.2 = c("Fibroblast_Day7"), test.use="MAST")
#Cluster 1
Endothelial_p5_DE_genes <- FindMarkers(all, ident.1 = c("Endothelial_Day0"), ident.2 = c("Endothelial_Day0p5"), test.use="MAST")
Endothelial_1_DE_genes <- FindMarkers(all, ident.1 = c("Endothelial_Day0"), ident.2 = c("Endothelial_Day1"), test.use="MAST")
Endothelial_4_DE_genes <- FindMarkers(all, ident.1 = c("Endothelial_Day0"), ident.2 = c("Endothelial_Day4"), test.use="MAST")
Endothelial_7_DE_genes <- FindMarkers(all, ident.1 = c("Endothelial_Day0"), ident.2 = c( "Endothelial_Day7"), test.use="MAST")
#Cluster 2
Monocyte_MacroPhage_DC_p5_DE_genes <- FindMarkers(all, ident.1 = c("Monocyte_MacroPhage_DC_Day0"), ident.2 = c("Monocyte_MacroPhage_DC_Day0p5"), test.use="MAST")
Monocyte_MacroPhage_DC_1_DE_genes <- FindMarkers(all, ident.1 = c("Monocyte_MacroPhage_DC_Day0"), ident.2 = c("Monocyte_MacroPhage_DC_Day1"), test.use="MAST")
Monocyte_MacroPhage_DC_4_DE_genes <- FindMarkers(all, ident.1 = c("Monocyte_MacroPhage_DC_Day0"), ident.2 = c("Monocyte_MacroPhage_DC_Day4"), test.use="MAST")
Monocyte_MacroPhage_DC_7_DE_genes <- FindMarkers(all, ident.1 = c("Monocyte_MacroPhage_DC_Day0"), ident.2 = c("Monocyte_MacroPhage_DC_Day7"), test.use="MAST")
#Cluster 3
B_Plasma_p5_DE_genes <- FindMarkers(all, ident.1 = c("B_Plasma_Day0"), ident.2 = c("B_Plasma_Day0p5"), test.use="MAST")
B_Plasma_1_DE_genes <- FindMarkers(all, ident.1 = c("B_Plasma_Day0"), ident.2 = c("B_Plasma_Day1"), test.use="MAST")
B_Plasma_4_DE_genes <- FindMarkers(all, ident.1 = c("B_Plasma_Day0"), ident.2 = c("B_Plasma_Day4"), test.use="MAST")
B_Plasma_7_DE_genes <- FindMarkers(all, ident.1 = c("B_Plasma_Day0"), ident.2 = c("B_Plasma_Day7"), test.use="MAST", min.cells.group = 2)


#Cluster 4
Epithelial_Mucin_p5_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial_Mucin_Day0"), ident.2 = c("Epithelial_Mucin_Day0p5"), test.use="MAST")
Epithelial_Mucin_1_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial_Mucin_Day0"), ident.2 = c("Epithelial_Mucin_Day1"), test.use="MAST")
Epithelial_Mucin_4_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial_Mucin_Day0"), ident.2 = c("Epithelial_Mucin_Day4"), test.use="MAST")
Epithelial_Mucin_7_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial_Mucin_Day0"), ident.2 = c("Epithelial_Mucin_Day7"), test.use="MAST")
#Cluster 5
# RBC_p5_DE_genes <- FindMarkers(all, ident.1 = c("RBC_Day0"), ident.2 = c("RBC_Day0p5"), test.use="MAST")
# RBC_1_DE_genes <- FindMarkers(all, ident.1 = c("RBC_Day0"), ident.2 = c("RBC_Day1"), test.use="MAST")
# RBC_4_DE_genes <- FindMarkers(all, ident.1 = c("RBC_Day0"), ident.2 = c("RBC_Day4"), test.use="MAST")
# RBC_7_DE_genes <- FindMarkers(all, ident.1 = c("RBC_Day0"), ident.2 = c("RBC_Day7"), test.use="MAST")
#Cluster 6
Lymphatic_Endothelial_Cell_p5_DE_genes <- FindMarkers(all, ident.1 = c("Lymphatic_Endothelial_Cell_Day0"), ident.2 = c("Lymphatic_Endothelial_Cell_Day0p5"), test.use="MAST")
Lymphatic_Endothelial_Cell_1_DE_genes <- FindMarkers(all, ident.1 = c("Lymphatic_Endothelial_Cell_Day0"), ident.2 = c("Lymphatic_Endothelial_Cell_Day1"), test.use="MAST")
Lymphatic_Endothelial_Cell_4_DE_genes <- FindMarkers(all, ident.1 = c("Lymphatic_Endothelial_Cell_Day0"), ident.2 = c("Lymphatic_Endothelial_Cell_Day4"), test.use="MAST")
Lymphatic_Endothelial_Cell_7_DE_genes <- FindMarkers(all, ident.1 = c("Lymphatic_Endothelial_Cell_Day0"), ident.2 = c("Lymphatic_Endothelial_Cell_Day7"), test.use="MAST")
#Cluster 7
T_p5_DE_genes <- FindMarkers(all, ident.1 = c("T _Day0"), ident.2 = c("T _Day0p5"), test.use="MAST")
T_1_DE_genes <- FindMarkers(all, ident.1 = c("T _Day0"), ident.2 = c("T _Day1"), test.use="MAST")
T_4_DE_genes <- FindMarkers(all, ident.1 = c("T _Day0"), ident.2 = c("T _Day4"), test.use="MAST")
T_7_DE_genes <- FindMarkers(all, ident.1 = c("T _Day0"), ident.2 = c("T _Day7"), test.use="MAST")
#Cluster 8
Epithelial_p5_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial_Day0"), ident.2 = c("Epithelial_Day0p5"), test.use="MAST")
Epithelial_1_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial_Day0"), ident.2 = c("Epithelial_Day1"), test.use="MAST")
Epithelial_4_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial_Day0"), ident.2 = c("Epithelial_Day4"), test.use="MAST")
Epithelial_7_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial_Day0"), ident.2 = c("Epithelial_Day7"), test.use="MAST")
#Cluster 9
Neutrophil_p5_DE_genes <- FindMarkers(all, ident.1 = c("Neutrophil_Day0"), ident.2 = c("Neutrophil_Day0p5"), test.use="MAST")
Neutrophil_1_DE_genes <- FindMarkers(all, ident.1 = c("Neutrophil_Day0"), ident.2 = c("Neutrophil_Day1"), test.use="MAST")
Neutrophil_4_DE_genes <- FindMarkers(all, ident.1 = c("Neutrophil_Day0"), ident.2 = c("Neutrophil_Day4"), test.use="MAST")
Neutrophil_7_DE_genes <- FindMarkers(all, ident.1 = c("Neutrophil_Day0"), ident.2 = c("Neutrophil_Day7"), test.use="MAST")
#Cluster 10
b_plasma_p5_DE_genes <- FindMarkers(all, ident.1 = c("b_plasma_Day0"), ident.2 = c("b_plasma_Day0p5"), test.use="MAST", min.cells.group = 2)
b_plasma_1_DE_genes <- FindMarkers(all, ident.1 = c("b_plasma_Day0"), ident.2 = c("b_plasma_Day1"), test.use="MAST")
b_plasma_4_DE_genes <- FindMarkers(all, ident.1 = c("b_plasma_Day0"), ident.2 = c("b_plasma_Day4"), test.use="MAST")
b_plasma_7_DE_genes <- FindMarkers(all, ident.1 = c("b_plasma_Day0"), ident.2 = c("b_plasma_Day7"), test.use="MAST")
#Cluster 11
Neural_p5_DE_genes <- FindMarkers(all, ident.1 = c("Neural_Day0"), ident.2 = c("Neural_Day0p5"), test.use="MAST")
Neural_1_DE_genes <- FindMarkers(all, ident.1 = c("Neural_Day0"), ident.2 = c("Neural_Day1"), test.use="MAST")
Neural_4_DE_genes <- FindMarkers(all, ident.1 = c("Neural_Day0"), ident.2 = c("Neural_Day4"), test.use="MAST")
Neural_7_DE_genes <- FindMarkers(all, ident.1 = c("Neural_Day0"), ident.2 = c("Neural_Day7"), test.use="MAST")
#Cluster 12
Smooth_Muscle_p5_DE_genes <- FindMarkers(all, ident.1 = c("Smooth_Muscle_Day0"), ident.2 = c("Smooth_Muscle_Day0p5"), test.use="MAST")
Smooth_Muscle_1_DE_genes <- FindMarkers(all, ident.1 = c("Smooth_Muscle_Day0"), ident.2 = c("Smooth_Muscle_Day1"), test.use="MAST")
Smooth_Muscle_4_DE_genes <- FindMarkers(all, ident.1 = c("Smooth_Muscle_Day0"), ident.2 = c("Smooth_Muscle_Day4"), test.use="MAST")
Smooth_Muscle_7_DE_genes <- FindMarkers(all, ident.1 = c("Smooth_Muscle_Day0"), ident.2 = c("Smooth_Muscle_Day7"), test.use="MAST")
#Cluster 13
epithelial_mucin_p5_DE_genes <- FindMarkers(all, ident.1 = c("epithelial_mucin_Day0"), ident.2 = c("epithelial_mucin_Day0p5"), test.use="MAST")
epithelial_mucin_1_DE_genes <- FindMarkers(all, ident.1 = c("epithelial_mucin_Day0"), ident.2 = c("epithelial_mucin_Day1"), test.use="MAST")
epithelial_mucin_4_DE_genes <- FindMarkers(all, ident.1 = c("epithelial_mucin_Day0"), ident.2 = c("epithelial_mucin_Day4"), test.use="MAST")
epithelial_mucin_7_DE_genes <- FindMarkers(all, ident.1 = c("epithelial_mucin_Day0"), ident.2 = c("epithelial_mucin_Day7"), test.use="MAST")
#Cluster 14
Skeletal_Muscle_p5_DE_genes <- FindMarkers(all, ident.1 = c("Skeletal_Muscle_Day0"), ident.2 = c("Skeletal_Muscle_Day0p5"), test.use="MAST")
Skeletal_Muscle_1_DE_genes <- FindMarkers(all, ident.1 = c("Skeletal_Muscle_Day0"), ident.2 = c("Skeletal_Muscle_Day1"), test.use="MAST")
Skeletal_Muscle_4_DE_genes <- FindMarkers(all, ident.1 = c("Skeletal_Muscle_Day0"), ident.2 = c("Skeletal_Muscle_Day4"), test.use="MAST")
Skeletal_Muscle_7_DE_genes <- FindMarkers(all, ident.1 = c("Skeletal_Muscle_Day0"), ident.2 = c("Skeletal_Muscle_Day7"), test.use="MAST", min.cells.group = 2)


DE_results<-grep("_DE_genes",names(.GlobalEnv),value=TRUE)
rm(DE_results_list)
DE_results_list<-do.call("list",mget(DE_results))

names(DE_results_list)<-sub("epithelial","Epi2",names(DE_results_list))
names(DE_results_list)<-sub("b","B2",names(DE_results_list))
names(DE_results_list)<-sub("thelial","",names(DE_results_list))
names(DE_results_list)<-sub("Monocyte_MacroPhage_DC","MMD",names(DE_results_list))

write.xlsx(DE_results_list, file = "~/20_015/outputs/All_celltypes_DE_results.xlsx",rowNames = TRUE)


