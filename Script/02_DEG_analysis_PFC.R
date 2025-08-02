###Differential expression analysis of PFC
## Andi Liu

## loading
library(Seurat)
library(Signac)
library(dplyr)
library(patchwork)
library(future)
library(stringr)
library(tidydr)
library(tidyverse)
library(viridis)
library(qs)
library(ggplot2)
library(readxl)

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

library(GenomicRanges)
library(GenomicFeatures)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)


## read data
object <- readRDS("cellbender_PFC_object_10.31.rds")
print(object)

Idents(object) <- "cluster_celltype"
table(Idents(object))

################################################################################
############################# modify the metadata ##############################
################################################################################
# pull out metadata
biospecimen_meta <- read.csv("./EOAD_biospecimen_metadata.csv")
meta <- object@meta.data

# add phenotype information
id <- match(meta$individual_ID,biospecimen_meta$Simple_ID)

meta$sex <- as.factor(biospecimen_meta[id,]$Sex)
meta$age <- biospecimen_meta[id,]$Age
meta$race <- as.factor(biospecimen_meta[id,]$Race)
meta$ethinicity <- as.factor(biospecimen_meta[id,]$Ethinicity)
meta$PMI <- biospecimen_meta[id,]$PMI.Hours.
meta$RIN <- biospecimen_meta[id,]$New_RIN
meta$batch <- as.factor(meta$batch)
#meta$batch <- biospecimen_meta[id,]$Batch

#meta

## assign back
object@meta.data <- meta

colnames(meta)

################################################################################
############ Find differentially expressed features between EOAD and NCI #######
################################################################################

##### filter genes
counts <- LayerData(object,assay = "PC",layer = "counts")# code in seurat 5.0+
genes.percent.expression <- rowMeans(counts>0 )*100   
genes.filter <- names(genes.percent.expression[genes.percent.expression>1])  #select genes expressed in at least 1% of cells

# using cell types with more than 50 cell in each batch and condition for now
cell.type_used <- c("Astrocyte","Excitatory","Inhibitory","Microglia","Oligodendrocyte","OPC")

# create empty dataframe to store data
condition_markers_cell.type <- data.frame()

# Performe cell type specific DE analysis using MAST method
# covariates to control: number of genes expressed, sex, age, PMI, RIN, batch
for (i in cell.type_used){
  print(i)
  # perform MAST DE analysis
  de.markers <- FindMarkers(object, ident.1 = "EOAD",ident.2 = "NCI", group.by = "diagnosis", # specify the conditions
                            subset.ident = i,assay = "PC",features =genes.filter,  # specify the cell type 
                            test.use = "MAST", # specify the method to use
                            min.pct = 0.25,
                            latent.vars = c("sex","age","batch")) # specify the covariates
  # filter the results based on the adjusted p value <= 0.05
  # de.markers <- de.markers[de.markers$p_val_adj <= 0.05,]
  # add cell type name
  de.markers$cell.type <- i
  # add gene to avoid duplication issue
  de.markers$gene <- rownames(de.markers)
  
  # add into the general data frame
  condition_markers_cell.type <- rbind(condition_markers_cell.type,de.markers)
  print(paste(i, "Done!"))
}

# check how many genes
table(condition_markers_cell.type$cell.type)

## adding direction information
condition_markers_cell.type$dir <- ifelse(condition_markers_cell.type$avg_log2FC>0,"pos","neg")
## check DEGs
table(condition_markers_cell.type$cell.type,condition_markers_cell.type$dir)
## saving results
write.csv(condition_markers_cell.type, file = "./Results/DEG/MAST_major.cell.type_PFC_7.19.csv")


## quitting
q("no")






