#PART 2: DOUBLET REMOVAL

library(Seurat)
library(scDblFinder)
set.seed(432)

##DAY0_UNL
LIG0_OLD <- readRDS("~/20_015/souped_DAY0_UNL.rds")
LIG0_OLDsce <- as.SingleCellExperiment(LIG0_OLD)
LIG0_OLDsce <- scDblFinder(LIG0_OLDsce)
LIG0_OLD <- as.Seurat(LIG0_OLDsce, counts = "counts", data = NULL)
Idents(LIG0_OLD) <- "scDblFinder.class"
LIG0_OLD <- subset(LIG0_OLD, idents = "singlet")
saveRDS(LIG0_OLD, file = "cleaned_DAY0_UNL.rds")

##DAY0
souped_LIG0 <- readRDS("~/20_015/souped_DAY0.rds")
LIG0sce <- as.SingleCellExperiment(souped_LIG0)
LIG0sce <- scDblFinder(LIG0sce)
LIG0 <- as.Seurat(LIG0sce, counts = "counts", data = NULL)
Idents(LIG0) <- "scDblFinder.class"
LIG0 <- subset(LIG0, idents = "singlet")
saveRDS(LIG0, file = "cleaned_DAY0.rds")

##DAY0.5
LIG5 <- readRDS("~/20_015/souped_DAY0.5.rds")
LIG5sce <- as.SingleCellExperiment(LIG5)
LIG5sce <- scDblFinder(LIG5sce)
LIG5 <- as.Seurat(LIG5sce, counts = "counts", data = NULL)
Idents(LIG5) <- "scDblFinder.class"
LIG5 <- subset(LIG5, idents = "singlet")
saveRDS(LIG5, file = "cleaned_DAY0.5.rds")

##DAY1_LIG
LIG1_OLD <- readRDS("~/20_015/souped_DAY1_LIG.rds")
LIG1_OLDsce <- as.SingleCellExperiment(LIG1_OLD)
LIG1_OLDsce <- scDblFinder(LIG1_OLDsce)
LIG1_OLD <- as.Seurat(LIG1_OLDsce, counts = "counts", data = NULL)
Idents(LIG1_OLD) <- "scDblFinder.class"
LIG1_OLD <- subset(LIG1_OLD, idents = "singlet")
saveRDS(LIG1_OLD, file = "cleaned_DAY1_LIG.rds")

##DAY1
LIG1 <- readRDS("~/20_015/souped_DAY1.rds")
LIG1sce <- as.SingleCellExperiment(LIG1)
LIG1sce <- scDblFinder(LIG1sce)
LIG1 <- as.Seurat(LIG1sce, counts = "counts", data = NULL)
Idents(LIG1) <- "scDblFinder.class"
LIG1 <- subset(LIG1, idents = "singlet")
saveRDS(LIG1, file = "cleaned_DAY1.rds")

##DAY4
LIG4 <- readRDS("~/20_015/souped_DAY4.rds")
LIG4sce <- as.SingleCellExperiment(LIG4)
LIG4sce <- scDblFinder(LIG4sce)
LIG4 <- as.Seurat(LIG4sce, counts = "counts", data = NULL)
Idents(LIG4) <- "scDblFinder.class"
LIG4 <- subset(LIG4, idents = "singlet")
saveRDS(LIG4, file = "cleaned_DAY4.rds")

##DAY7
LIG7 <- readRDS("~/20_015/souped_DAY7.rds")
LIG7sce <- as.SingleCellExperiment(LIG7)
LIG7sce <- scDblFinder(LIG7sce)
LIG7 <- as.Seurat(LIG7sce, counts = "counts", data = NULL)
Idents(LIG7) <- "scDblFinder.class"
LIG7 <- subset(LIG7, idents = "singlet")
saveRDS(LIG7, file = "cleaned_DAY7.rds")
