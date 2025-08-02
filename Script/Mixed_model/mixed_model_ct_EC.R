## DE analysis of multiome subjects using mixed linear model
# adapted from https://www.nature.com/articles/s41467-021-25960-2
# https://github.com/neurorestore/Libra/blob/main/R/mixedmodel_de.R
## Andi Liu
# 9/21/2023

library(Libra)
library(Seurat)
library(pbmcapply)

library(dplyr)# for dataframe processing
library(tidyr)#
library(tidydr)
library(magrittr) # for %<>%

library(glmmTMB)# negtive bionomial model
#library(lmerTest,lme4) # linear mixed model
#library(blme) #poisson model

# load previouse object
rm(list = ls())
gc()
object <- readRDS("./Data/cellbender_EC_object_1.9.rds")
object

################################################################################
################################ parameters ####################################
################################################################################
date() 
celltype <- "Astrocyte" #change to the cell type you want to analyze
min_features_pct <- 0.25
out_path <- "./Results/DEG/DE_mixed_EC_Ast.csv"

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
meta$PMI <- as.numeric(biospecimen_meta[id,]$PMI.Hours.)
meta$RIN <- biospecimen_meta[id,]$New_RIN
meta$batch <- as.factor(biospecimen_meta[id,]$Batch)

#meta
## assign back
object@meta.data <- meta
unique(meta$cluster_celltype)

################################################################################
####################################formula#####################################
################################################################################
fmla <- paste("GENE ~",paste("diagnosis", "sex","age","batch", sep = "+"))
fmla = as.formula(paste(fmla, "+ offset(log(total_counts)) +", "(1|replicate)"))

##null_fmla will be used using Likelihood ratio test
null_fmla = as.formula(gsub("diagnosis \\+ ", " ", deparse(fmla)))

fmla
null_fmla

################################################################################
###################################analysis#####################################
################################################################################
obj <- subset(object,subset = cluster_celltype == celltype)
print(celltype)
obj

################################################################################
################ running mixed model for DEG identification ####################
################################################################################
# still using counts from RNA assay
expr <- LayerData(obj,assay = "PC",layer = "counts")
#expr <- round(expr)
meta <- obj@meta.data
  
## check raw data
print(table(paste(meta$diagnosis, meta$individual_ID, sep=":")))
print(dim(expr))
print(dim(meta))
    
# check minimum features
genes.percent.expression <- rowMeans(expr>0)
#keep = rowSums(expr) >= min_features_pct*length(meta$major.cell.type.manual)
keep <-   names(genes.percent.expression[genes.percent.expression>min_features_pct])
expr = expr[keep,]

### check how many genes includes
print(dim(expr))
print(dim(meta))
genes = rownames(expr)
print(paste(length(genes),"genes will be tested."))

################################################################################
###################### ttttttttesssssssssstttttttttttting ######################
################################################################################
results <-pbmclapply(genes,mc.cores = 12,FUN = function(x){
    test_dat = data.frame(
      GENE = expr[x,],
      diagnosis = as.factor(meta$diagnosis),
      replicate = meta$individual_ID,
      sex = as.factor(meta$sex),
      age = meta$age,
      batch = as.factor(meta$batch),
      PMI = meta$PMI
    )
    gene = x
    # relevel diagnosis
    test_dat <- within(test_dat, diagnosis <- relevel(diagnosis, ref = "NCI"))
    ## adding offset
    test_dat %<>% mutate(total_counts = as.numeric(colSums(expr)+1))
    ## running test
    mod1 = glmmTMB(fmla, test_dat, family = nbinom1, REML = FALSE)
    mod2 = glmmTMB(null_fmla, test_dat, family = nbinom1, REML = FALSE)
    lrt = try(anova(mod1, mod2, refit = F))
    p_val = lrt$`Pr(>Chisq)`[2]
    test_statistic = lrt$Chisq[2]
    return(c(gene,p_val,test_statistic))
  })

## create final dataframe
final <- do.call(rbind, results)
colnames(final) <- c("genes", "p_val","test_statistic")
final <- as.data.frame(final)

# function for gettting logFC
run_logFC <- function(expr0,meta0){
  sc0 <-  CreateSeuratObject(expr0, meta = meta0) %>% NormalizeData()
  
  Idents(sc0) = relevel(as.factor(sc0$diagnosis),ref = "NCI")
  mat = LayerData(sc0, assay = "RNA",layer = 'data')
  levels = levels(as.factor(meta0$diagnosis))
  
  cells1 = WhichCells(sc0, idents = levels[1])
  cells2 = WhichCells(sc0, idents = levels[2])
  data1 = log(rowMeans(expm1(mat[, cells1, drop = F]) + 1),base = 2)
  data2 = log(rowMeans(expm1(mat[, cells2, drop = F]) + 1),base = 2)
  logFC = as.numeric(data1 - data2) # be careful about the direction, data1 (cases) - data2 (controls)
  
  return(logFC)
}

################################################################################
################ adding other variables and saving results #####################
################################################################################
final$avg_log2FC <- run_logFC(expr,meta)
final$p_val_adj = p.adjust(final$p_val, method = 'BH')
final$dir <- ifelse(final$avg_log2FC>0,"pos","neg")
final$cell.type <- celltype
final$DE_method <- "negbinom"
final$DE_type <- 'LRT'

write.csv(final, file = out_path)

q("no")