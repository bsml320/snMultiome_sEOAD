# loading
library(Seurat)
library(Signac)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(GenomicRanges)
library(GenomicFeatures)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

## read in related data
object <- readRDS("cellbender_HIP_object_10.31.rds")
object

DefaultAssay(object) <- "CTpeaks"

## normalization again in case
object <- RunTFIDF(object)
object <- FindTopFeatures(object, min.cutoff = 'q0')
object <- RunSVD(object)

object

################################################################################
######### Find differentially accessibility features between EOAD and NCI ######
################################################################################
meta <- object@meta.data
meta$age <- as.numeric(meta$age)
meta$sex <- as.character(meta$sex)
object@meta.data <- meta
# using cell types with more than 50 cell in each batch and condition for now
cell.type_used <- c("Astrocyte","Excitatory","Inhibitory","Microglia","Oligodendrocyte","OPC")

# create empty dataframe to store data
condition_markers_cell.type <- data.frame()

# Performe cell type specific DE analysis using MAST method
# covariates to control: number of genes expressed, sex, age, PMI, RIN, batch
for (i in cell.type_used){
  print(i)
  # perform MAST DE analysis
  da.markers <- suppressWarnings(FindMarkers(object, ident.1 = "EOAD",ident.2 = "NCI", group.by = "diagnosis",min.pct = 0.01, # specify the conditions
                                           subset.ident = i,assay = "CTpeaks",  # specify the cell type
                                           #features = unique(link_interested$peak),
                                           test.use = "LR", # specify the method to use
                                           latent.vars = c("batch","nCount_CTpeaks"))) # specify the covariates
  # filter the results based on the adjusted p value <= 0.05
  #da.markers <- da.markers[da.markers$p_val <= 0.05,]
  # add cell type name
  da.markers$cell.type <- i
  # add gene to avoid duplication issue
  da.markers$peaks <- rownames(da.markers)
  
  # add into the general data frame
  condition_markers_cell.type <- rbind(condition_markers_cell.type,da.markers)
  print(paste(i, "Done!"))
}

# check how many genes
table(condition_markers_cell.type$cell.type)
write.csv(condition_markers_cell.type, file = "./Results/DA/DA_major.cell.type_HIP.csv")

q("no")