### motif enrichment analysis of PFC/EC/HIP
## Andi Liu
# 12/13/2023
## loading
suppressPackageStartupMessages({
library(Seurat)
library(Signac)
library(dplyr)
library(patchwork)
library(future)
library(stringr)
library(tidydr)
library(tidyverse)
library(viridis)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# For motif analysis
library(JASPAR2020)
library(TFBSTools)

library('ComplexHeatmap')
library(circlize)
})

# Annotating the linked peaks
## loading packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# downloaded hg38 known gene annotation from UCSC: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/
# read into txdb format
#txdb <- makeTxDbFromGFF("./Data/hg38.knownGene.gtf.gz")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## read in related data
#object <- readRDS("cellbender_PFC_object_10.31.rds")
#object <- readRDS("cellbender_EC_object_10.31.rds")
object <- readRDS("cellbender_HIP_object_10.31.rds")
object


## download human CORE TFBS from JASPAR
fn <- file.path("./Data/JASPAR2024_hs_core_755.txt")
fn

pfm <- readJASPARMatrix(fn, matrixClass="PFM")
#names(pfm) <- str_split_fixed(names(pfm),"\\.",n = 3)[,3]
pfm

# build function to run motif enrichment by cell type
get_motif <- function(object,celltype,deg_direction){
    # get interest peaks
    peak_interest <- links[links$celltype == celltype & links$deg_dir == deg_direction & links$score >0.05,]$peak
    message(paste("Testing on",length(peak_interest),"peaks links to",celltype,deg_direction,"DEGs."))

    # find peaks open in selected cell type
    DefaultAssay(object) <- "CTpeaks"
    open.peaks <- AccessiblePeaks(object, idents = celltype)
    
    # match the overall GC content in the peak set
    meta.feature <- GetAssayData(object, assay = "CTpeaks", slot = "meta.features")
    peaks.matched <- MatchRegionStats(
        meta.feature = meta.feature[open.peaks, ],
        query.feature = meta.feature[peak_interest, ],
        n = 50000)

    ## test enrichment
    enriched.motifs <- FindMotifs(
        object = object,
        features = peak_interest,
        background=peaks.matched)
    
    enriched.motifs$celltype <- celltype
    enriched.motifs$deg_dir <- deg_direction
    
    # getting significant results
    sig.enriched.motifs <- enriched.motifs[enriched.motifs$p.adjust < 0.05,]
    
    # results
    return(sig.enriched.motifs)
}

# get celltype specific peaks
DefaultAssay(object) <- "CTpeaks"
# add motif information
object <- AddMotifs(object = object,genome = BSgenome.Hsapiens.UCSC.hg38,pfm = pfm)

## Get mortif information for finding TF - peaks - differential correlation - DEGs
motif.all <- GetMotifData(
    object = object, assay = "CTpeaks", slot = "data"
  )
# motif.all[1,]
motif.names <- GetMotifData(
    object = object, assay = "CTpeaks", slot = "motif.names"
  )
length(motif.names)
#names(motif.names$V1) <- str_split_fixed(names(motif.names),"\\.",n = 3)[,3]
# generate dataframe of motif name and presents.
motif.names <- as.data.frame(as.matrix(motif.names))
motif.names$V1 <- str_split_fixed(motif.names$V1,"\\.",n = 3)[,3]
motif.names$V2 <- str_split_fixed(motif.names$V1,"::",n = 2)[,1]

head(motif.names)
#unlist(motif.names)
#rownames(motif.names[motif.names$V1 == "ASCL1",])
table(unique(motif.names$V1) %in% rownames(object@assays$PC))
table(unique(motif.names$V2) %in% rownames(object@assays$PC))

length(unique(motif.names$V2))


## running analysis on PFC region object
# load data
object

## top upregulated genes
degs <- read.csv("./Results/DEG/Overlap_mast_mixed_HIP.csv",row.names = 1)
message("How many DEGs in each cell type in HIP.")
dim(degs)
table(degs$celltype,degs$dir)
#degs

## read in the link data
links<- read.csv("./Results/LINK/HIP_linkpeaks_all_annotated_1.23.csv")
links$new_comb <- paste(links$gene,links$deg_dir,links$celltype,sep = "_")
links <- links[links$width > 3000,]
#links <- links[links$p.adj < 0.05,]
#links <- links[links$score > 0.1 & links$p.adj < 0.05,]
#
message("How many linked peaks in each cell type in PFC.")
dim(links)
table(links$celltype,links$deg_dir)
#
message("How many linked DEGs in each cell type in PFC.")
degs <- degs[degs$comb %in% unique(links$new_comb),]
table(degs$celltype,degs$dir)
length(degs$comb)
#degs
links <- links[links$new_comb %in% degs$comb,]

### running for each cell type and direction
ct <- c("Astrocyte","Excitatory","Inhibitory","Microglia","Oligodendrocyte","OPC")
direction <- c("pos","neg")
Idents(object) <- "cluster_celltype"
final <- data.frame()

#i = "Excitatory"
#j = "neg" 
for (i in ct){
    for(j in direction){
        celltype <- i
        print(celltype)
        deg_direction = j
        print(dir)
        
        res1 <- get_motif(object=object, celltype=celltype, deg_direction = deg_direction)
        print(length(res1$motif.name))
        final <- rbind(final,res1)
    }
}

#### generate heatmap for motif enriched results_PFC
enriched.motif <- final
table(enriched.motif$celltype,enriched.motif$deg_dir)

enriched.motif$motif.name <- str_split_fixed(enriched.motif$motif.name,"\\.",n = 3)[,3]
enriched.motif$motif.name <- str_split_fixed(enriched.motif$motif.name,"::",n = 1)[,1]

head(enriched.motif)

## only keeping those TFs expressed in the specific cell type. 
expr.res <- data.frame()
ct = c("Astrocyte","Excitatory","Inhibitory","Microglia","Oligodendrocyte","OPC")

for (i in ct){
    print(i)

    # calculating co-expression matrix for cell type of interest
    obj <- subset(object, subset = cluster_celltype == i)
    
    # only focus on expressed protein coding genes in the cell type of interes
    counts <- LayerData(obj,assay = "PC",layer = "counts")# code in seurat 5.0+
    genes.percent.expression <- rowMeans(counts>0 )*100   
    genes_selected <- names(genes.percent.expression[genes.percent.expression>25])
    message(paste(length(genes_selected)," genes expressed in >25% cells."))

    motif.res <- enriched.motif[enriched.motif$celltype == i,]
    motif.res <- motif.res[motif.res$motif.name %in% genes_selected,]
    message(length(motif.res$motif))

    expr.res <-rbind(expr.res,motif.res)

}

table(expr.res$celltype,expr.res$deg_dir)
## saving results
write.csv(expr.res, file = "./Results/LINK/HIP_enriched_JASPAR_motif_expr_1.23.csv")
## saving files
saveRDS(object, file = "cellbender_HIP_object_10.31.rds")