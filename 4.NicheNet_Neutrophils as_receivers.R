library(nichenetr)
library(tidyverse)
library(Seurat)
library(future)
library(ggplot2)
plan("multisession", workers = 4)
options(future.globals.maxSize = 500 * 1024^2)

setwd("~/20_015/nichenet_output3/")
set.seed(42)


##Load NicheNet Database Method 1: automatic read-in from internet
 # organism <- "mouse"
 # 
 # if(organism == "human"){
 #   lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
 #   ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
 #   weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
 # } else if(organism == "mouse"){
 #   lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
 #   ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
 #   weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
 #   
 # }


##Load NicheNet Database Method 2: -manual: download & read-in from local computer-
ligand_target_matrix = readRDS("~/20_015/nichenet_data/nichenetR/ligand_target_matrix_nsga2r_final_mouse.rds")
lr_network = readRDS("~/20_015/nichenet_data/nichenetR/lr_network_mouse_21122021.rds")
weighted_networks = readRDS("~/20_015/nichenet_data/nichenetR/weighted_networks_nsga2r_final_mouse.rds")



##Read saved RDS from 3.Analysis (FindMarkers Analysis)
seuratObj <- readRDS("~/20_015/outputs/res0.1/all_annotated_without_RBC_Skeletal.RDS")

table(Idents(seuratObj))

seuratObj$celltype<-Idents(seuratObj)
seuratObj@meta.data$celltype %>% table()
seuratObj@meta.data$tissue %>% table()
table(seuratObj$celltype,seuratObj$tissue)

Idents(seuratObj)<-"celltype"
table(Idents(seuratObj))



##Set NicheNet Parameters: - manual: this condition compares Fibroblasts_Day0 to Fibroblast_Day0p5 against all other celltypes - 
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver =  c("Neutrophil"),
  condition_colname = "tx_time",
  condition_oi = c("Day0p5","Day1"),
  condition_reference = "Day0", 
  sender = c("Endothelial", "B", "Mononuclear_Phagocyte", "Epithelial", "Fibroblast","T", "Neural", "Smooth_Muscle"),
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks,
  top_n_ligands = 20,
  top_n_targets = 50,
  lfc_cutoff = 0.25)



##DE Features, Ligands, Receptor Figures:
levels(seuratObj)
nichenet_output$ligand_expression_dotplot
nichenet_output$ligand_activities
nichenet_output$top_ligands
nichenet_output$top_targets


#Ligand vs DE Features Regulatory Potential Heatmap: 
  pdf(file="~/20_015/nichenet_output3/DE_genes_with_ligands.pdf", width = 7, height = 4.5)
  nichenet_output$ligand_target_heatmap +
   scale_fill_gradient2(low = "whitesmoke",  high = "blue", breaks = c(0,0.0045,0.009)) +
   xlab ("DE genes between Gingiva and ligature fibroblasts")  +
   ylab("Prioritized cell ligands")
  dev.off()


#Unused Figures:
  nichenet_output$ligand_target_matrix %>% .[1:8,1:6]
  nichenet_output$ligand_target_df # weight column = regulatory potential
  nichenet_output$top_targets

  #pdf(file="Doplot_ligands_by_celltype_condition.pdf", width = 10, height = 3)
  #DotPlot(seuratObj %>% subset(idents = c("endothelial","epithelial","smooth_muscle","lymphatic_endothelial")), features = nichenet_output$top_ligands %>% rev(), split.by = "tissue", cols = c("blue","darkorange3")) + RotatedAxis()
  #dev.off()

  #VlnPlot(seuratObj %>% subset(idents = c("endothelial","epithelial","smooth_muscle","lymphatic_endothelial")), features = c("VCAN"), split.by = "tissue", pt.size = 0, combine = TRUE)

  
  
#Ligand Activity Summary Figures:
pdf(file="~/20_015/nichenet_output3/summary_ligand_activity_target_heatmap_fibroblast_target.pdf", width = 14, height = 8)
nichenet_output$ligand_activity_target_heatmap
dev.off()

pdf(file="~/20_015/nichenet_output3/summary_ligand_activity_target_heatmap_fibroblast_target_tall.pdf", width = 14, height = 9)
nichenet_output$ligand_activity_target_heatmap
dev.off()


#Ligand vs Receptor Heatmap: 
pdf(file="~/20_015/nichenet_output3/fibroblast_receptor.pdf", width = 5.75, height = 6)
nichenet_output$ligand_receptor_heatmap
dev.off()

nichenet_output$ligand_receptor_matrix %>% .[1:10,1:6]
nichenet_output$ligand_receptor_df
nichenet_output$top_receptors




## Ligand-Receptor Circos
# Assign ligands to sender cells
table(Idents(seuratObj))

avg_expression_ligands = AggregateExpression(seuratObj, features = nichenet_output$top_ligands, group.by = "celltype")

sender_ligand_assignment = avg_expression_ligands$SCT %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()

sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment) 

all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = nichenet_output$top_ligands %>% setdiff(unique_ligands)
Fibroblasts_specific_ligands = sender_ligand_assignment$Fibroblastibroblasts %>% names() %>% setdiff(general_ligands)
Endothelial_specific_ligands = sender_ligand_assignment$Endothelial %>% names() %>% setdiff(general_ligands)
B_specific_ligands = sender_ligand_assignment$B %>% names() %>% setdiff(general_ligands)
`Mononuclear-Phagocyte_specific_ligands` = sender_ligand_assignment$`Mononuclear-Phagocyte` %>% names() %>% setdiff(general_ligands)
Epithelial_specific_ligands = sender_ligand_assignment$Epithelial %>% names() %>% setdiff(general_ligands)
Neural_specific_ligands = sender_ligand_assignment$Neural %>% names() %>% setdiff(general_ligands)
T_specific_ligands = sender_ligand_assignment$T %>% names() %>% setdiff(general_ligands)


ligand_type_indication_df = tibble(
  ligand_type = c(rep("Fibroblast_specific", times = Fibroblasts_specific_ligands %>% length()),
                  rep("Endothelial_specific", times = Endothelial_specific_ligands %>% length()),
                  rep("B_specific", times = B_specific_ligands %>% length()),
                  rep("Mononuclear-Phagocyte_specific", times = `Mononuclear-Phagocyte_specific_ligands` %>% length()),
                  rep("Epithelial_specific", times = Epithelial_specific_ligands %>% length()),
                  rep("Neural_specific", times = Neural_specific_ligands %>% length()),
                  rep("T_specific", times = T_specific_ligands %>% length()),
                  rep("general", times = general_ligands %>% length())),
  
  ligand = c(Fibroblasts_specific_ligands,
             Endothelial_specific_ligands,
             B_specific_ligands,
             `Mononuclear-Phagocyte_specific_ligands`,
             Epithelial_specific_ligands,
             Neural_specific_ligands,
             T_specific_ligands,
             general_ligands))



#ligand_type_indication_df <- assign_ligands_to_celltype(seuratObj, nichenet_output$top_ligands, celltype_col = "celltypes", slot = "data")
ligand_type_indication_df %>% head()
ligand_type_indication_df$ligand_type %>% table()

# Define the ligand-target links of interest
head(nichenet_output$ligand_target_df)

active_ligand_target_links_df <- nichenet_output$ligand_target_df
active_ligand_target_links_df$target_type <- "Gingiva_vs_Perio-Neutrophil_DE" # needed for joining tables
circos_links <- get_ligand_target_links_oi(ligand_type_indication_df,
                                           active_ligand_target_links_df,
                                           cutoff = 0)


head(circos_links)

#Ligand Target Circos
ligand_colors <- c("general" = "green3",
                   "Endothelial_specific" = "purple") 

target_colors <- c("Gingiva_vs_Perio-Neutrophil_DE" = "steelblue2") 

vis_circos_obj <- prepare_circos_visualization(circos_links,
                                               ligand_colors = ligand_colors,
                                               target_colors = target_colors,
                                               celltype_order = c("Endothelial_specific",
                                                                  "general"))


pdf(file = "~/20_015/nichenet_output3/Circos_Perio_Neutrophils_Ligands-Targets.pdf", width = 6, height = 6)
make_circos_plot(vis_circos_obj, transparency = TRUE,  args.circos.text = list(cex = 0.5)) 
dev.off()


## If a legend is needed
par(bg = "transparent")

# Default celltype order
celltype_order <- unique(circos_links$ligand_type) %>% sort() %>% .[. != "General"] %>% c(., "General")

# Create legend
circos_legend <- ComplexHeatmap::Legend(
  labels = celltype_order,
  background = ligand_colors[celltype_order],
  type = "point",
  grid_height = unit(3, "mm"),
  grid_width = unit(3, "mm"),
  labels_gp = grid::gpar(fontsize = 8)
)

circos_legend_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(circos_legend))

make_circos_plot(vis_circos_obj, transparency = TRUE, args.circos.text = list(cex = 0.5))
p_circos_no_legend <- recordPlot()

cowplot::plot_grid(p_circos_no_legend, circos_legend_grob, rel_widths = c(1, 0.2))






#Ligand-Receptor Circos: 
lr_network_top_df <- nichenet_output$ligand_receptor_df %>%
  mutate(target_type = "Neutrophil_receptors") %>%
  rename(target=receptor) %>%
  inner_join(ligand_type_indication_df)

receptor_colors <- c("Neutrophil_receptors" = "gray")



table(circos_links$ligand_type)

vis_circos_receptor_obj <- prepare_circos_visualization(lr_network_top_df,
                                                        ligand_colors = c("Endothelial_specific" = "cyan2","Neural_specific" = "purple", "Mononuclear-Phagocyte_specific" = "black", "general" = "green3"),
                                                        target_colors = receptor_colors,
                                                        celltype_order = c("Endothelial_specific",
                                                                           "Neural_specific",
                                                                           "general"))

pdf(file = "~/20_015/nichenet_output3/Circos_Perio_Neutrophil_Ligands-Receptors.pdf", width = 8, height = 8)
make_circos_plot(vis_circos_receptor_obj,
                 transparency = TRUE, 
                 link.visible = TRUE,args.circos.text = list(cex = 0.8)) 
dev.off()




