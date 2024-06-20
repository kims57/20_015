### PART 3: Analysis

library(Seurat)
library(ggplot2)
library(clustree)
library(patchwork)
library(openxlsx)

setwd("~/20_015/data/")


## Load Datasets & Group by Tx Time & Percent Mitochondrial Filtering -manual: select percent mt threshold-
Day0 <- readRDS("~/20_015/data/cleaned_Day0.RDS")
Day0$tx_time <- "Day0"
Day0 <- PercentageFeatureSet(Day0, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(Day0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Day0 <- subset(Day0, subset = nFeature_RNA >= 200 & percent.mt < 10)

Day0p5 <- readRDS("~/20_015/data/cleaned_Day0.5.RDS")
Day0p5$tx_time <- "Day0p5"
Day0p5 <- PercentageFeatureSet(Day0p5, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(Day0p5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Day0p5 <- subset(Day0p5, subset = nFeature_RNA >= 200 & percent.mt < 8)

Day1 <- readRDS("~/20_015/data/cleaned_Day1.RDS")
Day1$tx_time <- "Day1"
Day1 <- PercentageFeatureSet(Day1, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(Day1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Day1 <- subset(Day1, subset = nFeature_RNA >= 200 & percent.mt < 7)

Day4 <- readRDS("~/20_015/data/cleaned_Day4.RDS")
Day4$tx_time <- "Day4"
Day4 <- PercentageFeatureSet(Day4, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(Day4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Day4 <- subset(Day4, subset = nFeature_RNA >= 200 & percent.mt < 5)

Day7 <- readRDS("~/20_015/data/cleaned_Day7.RDS")
Day7$tx_time <- "Day7"
Day7 <- PercentageFeatureSet(Day7, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(Day7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Day7 <- subset(Day7, subset = nFeature_RNA >= 200 & percent.mt < 5)


## Straight Merge 
all <- merge(x= Day0 , y = c(Day0p5, Day1, Day4, Day7), add.cell.ids = c("Day0","Day0p5","Day1","Day4","Day7"), project = "all")
saveRDS(all, file = "~/20_015/outputs/straight_merged_all.rds")
DefaultAssay(all)


## SC Transformation
all <- SCTransform(all, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = TRUE)


## Linear Dimensional Reduction (PCA)
all <- RunPCA(all, features = VariableFeatures(object = all))
DimPlot(all, reduction = "pca", shuffle = TRUE)
ElbowPlot(all, ndims=50)


##Optimal Clustering Analysis [File Output] -manual:select clustering resolution-
all <- FindNeighbors(all, dims = 1:50)
all <- FindClusters(all, resolution = c(0,0.1,0.2,0.3,0.4,0.5))
all <- RunUMAP(all, dims = 1:50)

pdf('~/20_015/outputs/all_clustree_seurat.pdf', width = 30, height = 20)
clustree(all, node_size_range=c(10,20), node_text_size = 8)
dev.off()

#Final Cluster Resolution Selection
Idents(all)<-"SCT_snn_res.0.1" 
all <- FindClusters(all, resolution=0.1, cluster.name = "seurat_clusters") 


## Differentially Expressed Features [File Output]
all_DE_markers <- FindAllMarkers(all, logfc.threshold = 0.5,min.pct = 0.25)
write.csv(all_DE_markers,"~/20_015/outputs/all_DE_markers.csv", row.names = FALSE)




##Cell Annotation Here -programatic or manual(EnrichR, Immugen, Panglao)-




##Assigning Cell Identity to Clusters
new.cluster.ids <- c("Fibroblast", "Endothelial", "B", "Mononuclear_Phagocyte", "Epithelial", "RBC",  "Endothelial", "Epithelial", "Neutrophil","T", "Neural", "Epithelial3", "Smooth_Muscle", "Skeletal_Muscle")
names(new.cluster.ids) <- levels(all)
all <- RenameIdents(all, new.cluster.ids)

all$celltype <- Idents(all)



##Save Original Project Prior to Edit [File Output]
saveRDS(all, file = "~/20_015/outputs/backup_all_RBC.RDS")


## Remove RBC Cluster
all <- subset(all, idents = c("RBC", "Skeletal_Muscle"), invert=TRUE)


# Cell Selector Debris Removal -manual:select cells to remove-
DimPlot(all, reduction="umap", group.by = "celltype", shuffle = TRUE)
p <- DimPlot(all, reduction="umap", group.by = "celltype", shuffle = TRUE)
selected_cells<-CellSelector(p)

all <- subset(all, cells =  selected_cells, invert=TRUE)



##Figures [File Output]
pdf(file="~/20_015/outputs/res0.1/Fig0_UMAP_Celltype_label.pdf")
DimPlot(all, reduction="umap", group.by = "celltype", shuffle = TRUE, label=TRUE) + xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("Cell Type")
dev.off()

pdf(file="~/20_015/outputs/res0.1/Fig1_UMAP_Celltype.pdf")
DimPlot(all, reduction="umap", group.by = "celltype", shuffle = TRUE) + xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("Cell Type")
dev.off()

pdf(file="~/20_015/outputs/res0.1/Fig2_UMAP_Txt_Time.pdf", width = 10, height = 8)
DimPlot(all, reduction = "umap", group.by = "tx_time", shuffle = TRUE)+ xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("Txt Time")
dev.off()

pdf(file="~/20_015/outputs/res0.1/Fig3_UMAP_Txt_Time_Split.pdf", width = 10, height = 8)
DimPlot(all, reduction = "umap", split.by = "tx_time", shuffle = TRUE)+ xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("Split Txt Time")
dev.off()


##Cell Proportion Figure [File Output]: 
table(Idents(all))
#prop.table(table, 1)
#prop.table(mytable, 2)

pdf(file='~/20_015/outputs/res0.1/Fig4_proportions.pdf')
mytable <- table(all$celltype, all$tx_time)
mydf <- as.data.frame(prop.table(mytable, 2))
ggplot(mydf, aes(x=Var2, y = Freq, fill=Var1)) + 
  geom_bar(colour="black", stat="identity") +
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




##DE Feature DotPlot Figure: [File Output] -manual: use DE markers to select features-
  #DotPlot(all, features = c("Col1a1", "Col1a2", "Col5a2", "Cd68", "Lyz2", "Fcer1g", "S100a8", "S100a9", "Il1b", "Igkc", "H2-Aa", "Ighm", "Trbc2", "Cd3g"), group.by = "celltype") + RotatedAxis()
  #DotPlot(all, features = c(""), group.by = "celltype") + RotatedAxis()



## FindMarkers: DE analysis btwn celltype_day0 vs celltype_dayx [Output File]
all$celltype.timepoint <- paste(all$celltype,all$tx_time, sep = "_")
Idents(all)<-"celltype.timepoint"

table(all$celltype, all$tx_time)

#FindMarkers Method 1: For Loop
# sample <- c("Fibroblast", "Endothelial", "B", "Mononuclear_Phagocyte", "Epithelial", "RBC",  "Endothelial", "Epithelail", "Neutrophil","T", "Neural", "Epithelial3", "Smooth_Muscle", "Skeletal_Muscle")
# for(x in sample) {
#   p5_DE_genes <- FindMarkers(all, ident.1 = c(paste0(x, "_Day0")), ident.2 = c(paste0(x, "_Day0p5")), test.use="MAST")
#   write.csv(p5_DE_genes, paste0(x,"_p5_DE_markers.csv"), row.names = FALSE)  # 
#   
#   one_DE_genes <- FindMarkers(all, ident.1 = c(paste0(x, "_Day0")), ident.2 = c( paste0(x, "_Day1")), test.use="MAST")
#   write.csv(one_DE_genes, paste0(x,"_1_DE_markers.csv"), row.names = FALSE)
#   
#   four_DE_genes <- FindMarkers(all, ident.1 = c(paste0(x, "_Day0")), ident.2 = c( paste0(x, "_Day4")), test.use="MAST")
#   write.csv(four_DE_genes, paste0(x,"_4_DE_markers.csv"), row.names = FALSE)
#   
#   seven_DE_genes <- FindMarkers(all, ident.1 = c(paste0(x, "_Day0")), ident.2 =  c(paste0(x, "_Day7")), test.use="MAST")
#   write.csv(seven_DE_genes, paste0(x,"_7_DE_markers.csv"), row.names = FALSE)
#   
#   rm(p5_DE_genes, one_DE_genes, four_DE_genes, seven_DE_genes)
#   gc()
# }


#FindMarkers Method 2: 
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
  B_p5_DE_genes <- FindMarkers(all, ident.1 = c("B_Day0"), ident.2 = c("B_Day0p5"), test.use="MAST")
  B_1_DE_genes <- FindMarkers(all, ident.1 = c("B_Day0"), ident.2 = c("B_Day1"), test.use="MAST")
  B_DE_genes <- FindMarkers(all, ident.1 = c("B_Day0"), ident.2 = c("B_Day4"), test.use="MAST")
  B_DE_genes <- FindMarkers(all, ident.1 = c("B_Day0"), ident.2 = c("B_Day7"), test.use="MAST")
  #Cluster 3
  Mononuclear_Phagocyte_p5_DE_genes <- FindMarkers(all, ident.1 = c("Mononuclear_Phagocyte_Day0"), ident.2 = c("Mononuclear_Phagocyte_Day0p5"), test.use="MAST")
  Mononuclear_Phagocyte_1_DE_genes <- FindMarkers(all, ident.1 = c("Mononuclear_Phagocyte_Day0"), ident.2 = c("Mononuclear_Phagocyte_Day1"), test.use="MAST")
  Mononuclear_Phagocyte_4_DE_genes <- FindMarkers(all, ident.1 = c("Mononuclear_Phagocyte_Day0"), ident.2 = c("Mononuclear_Phagocyte_Day4"), test.use="MAST")
  Mononuclear_Phagocyte_7_DE_genes <- FindMarkers(all, ident.1 = c("Mononuclear_Phagocyte_Day0"), ident.2 = c("Mononuclear_Phagocyte_Day7"), test.use="MAST")
  #Cluster 4
  Epithelial_p5_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial_Day0"), ident.2 = c("Epithelial_Day0p5"), test.use="MAST")
  Epithelial_1_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial_Day0"), ident.2 = c("Epithelial_Day1"), test.use="MAST")
  Epithelial_4_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial_Day0"), ident.2 = c("Epithelial_Day4"), test.use="MAST")
  Epithelial_7_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial_Day0"), ident.2 = c("Epithelial_Day7"), test.use="MAST")
  #Cluster 6 RBC
  #Cluster 8
  Neutrophil_p5_DE_genes <- FindMarkers(all, ident.1 = c("Neutrophil_Day0"), ident.2 = c("Neutrophil_Day0p5"), test.use="MAST")
  Neutrophil_1_DE_genes <- FindMarkers(all, ident.1 = c("Neutrophil_Day0"), ident.2 = c("Neutrophil_Day1"), test.use="MAST")
  Neutrophil_4_DE_genes <- FindMarkers(all, ident.1 = c("Neutrophil_Day0"), ident.2 = c("Neutrophil_Day4"), test.use="MAST")
  Neutrophil_7_DE_genes <- FindMarkers(all, ident.1 = c("Neutrophil_Day0"), ident.2 = c("Neutrophil_Day7"), test.use="MAST")
  #Cluster 9
  T_p5_DE_genes <- FindMarkers(all, ident.1 = c("T_Day0"), ident.2 = c("T_Day0p5"), test.use="MAST")
  T_1_DE_genes <- FindMarkers(all, ident.1 = c("T_Day0"), ident.2 = c("T_Day1"), test.use="MAST")
  T_4_DE_genes <- FindMarkers(all, ident.1 = c("T_Day0"), ident.2 = c("T_Day4"), test.use="MAST")
  T_7_DE_genes <- FindMarkers(all, ident.1 = c("T_Day0"), ident.2 = c("T_Day7"), test.use="MAST")
  #Cluster 10
  Neural_p5_DE_genes <- FindMarkers(all, ident.1 = c("Neural_Day0"), ident.2 = c("Neural_Day0p5"), test.use="MAST", min.cells.group = 2)
  Neural_1_DE_genes <- FindMarkers(all, ident.1 = c("Neural_Day0"), ident.2 = c("Neural_Day1"), test.use="MAST")
  Neural_4_DE_genes <- FindMarkers(all, ident.1 = c("Neural_Day0"), ident.2 = c("Neural_Day4"), test.use="MAST")
  Neural_7_DE_genes <- FindMarkers(all, ident.1 = c("Neural_Day0"), ident.2 = c("Neural_Day7"), test.use="MAST")
  #Cluster 11
  Epithelial3_p5_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial3_Day0"), ident.2 = c("Epithelial3_Day0p5"), test.use="MAST", min.cells.group = 1)
  Epithelial3_1_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial3_Day0"), ident.2 = c("Epithelial3_Day1"), test.use="MAST", min.cells.group = 1)
  #Epithelial3_4_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial3_Day0"), ident.2 = c("Epithelial3_Day4"), test.use="MAST", min.cells.group = 1) #not enough cells
  #Epithelial3_7_DE_genes <- FindMarkers(all, ident.1 = c("Epithelial3_Day0"), ident.2 = c("Epithelial3_Day7"), test.use="MAST", min.cells.group = 1) #not enough cells
  #Cluster 12
  Smooth_Muscle_p5_DE_genes <- FindMarkers(all, ident.1 = c("Smooth_Muscle_Day0"), ident.2 = c("Smooth_Muscle_Day0p5"), test.use="MAST")
  Smooth_Muscle_1_DE_genes <- FindMarkers(all, ident.1 = c("Smooth_Muscle_Day0"), ident.2 = c("Smooth_Muscle_Day1"), test.use="MAST")
  Smooth_Muscle_4_DE_genes <- FindMarkers(all, ident.1 = c("Smooth_Muscle_Day0"), ident.2 = c("Smooth_Muscle_Day4"), test.use="MAST")
  Smooth_Muscle_7_DE_genes <- FindMarkers(all, ident.1 = c("Smooth_Muscle_Day0"), ident.2 = c("Smooth_Muscle_Day7"), test.use="MAST")
  #Cluster 13: Skeletal Muscle

#save complied FindMarkers comparisons excel file [File Output]
DE_results<-grep("_DE_genes",names(.GlobalEnv),value=TRUE)
rm(DE_results_list)
DE_results_list<-do.call("list",mget(DE_results))

  #shorten title due to character limit: -manual: replace "thelial" w/ "" - 
  names(DE_results_list)<-sub("thelail","",names(DE_results_list))
  names(DE_results_list)<-sub("thelial","",names(DE_results_list))
  names(DE_results_list)<-sub("_Muscle","",names(DE_results_list))
  names(DE_results_list)<-sub("lear_Phagocyte","",names(DE_results_list))
  
write.xlsx(DE_results_list, file = "~/20_015/outputs/res0.1/FindMarkers_all_DE_results.xlsx",rowNames = TRUE)




##Save Analysis: 
saveRDS(all, file = "~/20_015/outputs/res0.1/all_annotated_without_RBC_Skeletal.RDS")




