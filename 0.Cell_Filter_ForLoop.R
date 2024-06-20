### PART 1: scRNA Filtering For Loop

library(Seurat)
library(SoupX)
library(stringr)




#Pathway Setup
setwd("~/20_015/data/")
samples <-c("Day0","Day0_UNL","Day0.5","Day1","Day1_LIG","Day4", "Day7")


# Cell-Cycle Scoring & Regression Gene List
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

#Change Human Gene to Mouse Gene Nomenclature
g2m.genes <- str_to_title(g2m.genes) 
s.genes <- str_to_title(s.genes) 




for (element in samples) {
 
  # SoupX Ambient RNA Removal
  tod <- Read10X_h5(filename =  paste0("~/20_015/",element,"/raw_feature_bc_matrix.h5"))
  toc <- Read10X_h5(filename = paste0("~/20_015/",element,"/filtered_feature_bc_matrix.h5"))
  sc = SoupChannel(tod, toc)
  srat = CreateSeuratObject(sc$toc)
  srat <- NormalizeData(srat)
  srat <- FindVariableFeatures(srat, nfeatures = 4000)
  srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat <- ScaleData(srat, vars.to.regress = c("S.Score", "G2M.Score"), verbose = TRUE)
  srat <- RunPCA(srat, npcs = 50, verbose = TRUE)
  srat <- FindNeighbors(srat, reduction = "pca", dims = 1:50)
  srat <- FindClusters(srat, resolution = 1)
  sc <- setClusters(sc,srat@active.ident)
  sc <- autoEstCont(sc)
  sample.deSouped <- adjustCounts(sc)
  sample <- CreateSeuratObject(sample.deSouped, project = element)

  
  # doublet id and removal
  library(scDblFinder)
  set.seed(432)
  sample_sce <- as.SingleCellExperiment(sample)
  sample_sce <- scDblFinder(sample_sce)
  sample_Dbl <- as.Seurat(sample_sce, counts = "counts", data = NULL)
  Idents(sample_Dbl) <- "scDblFinder.class"
  sample_clean <- subset(sample_Dbl, idents = "singlet")
  saveRDS(sample_clean, file = paste0("cleaned_", element,".RDS"))
  rm(sample,sample_sce,sample_Dbl, sample_clean)
  gc()
}