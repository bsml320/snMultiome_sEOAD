# runnning Peak to gene links analysis on each brain regions
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
# library for plotting
library(ggplot2)
library(dplyr)
library(viridis)
library(future)

# Annotating the linked peaks
## loading packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
options(future.globals.maxSize = 800 * 1024^3)
# downloaded hg38 known gene annotation from UCSC: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/
# read into txdb format
#txdb <- makeTxDbFromGFF("./Data/hg38.knownGene.gtf.gz")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## read in related data
object <- readRDS("cellbender_EC_object_10.31.rds")
object

## read in DEG list
DEG_list <- read.csv("./Results/DEG/Overlap_mast_mixed_EC.csv")
table(DEG_list$celltype,DEG_list$dir)
table(DEG_list$dir)

## build find link pipeline
get_link <- function(object, celltype, dir, DEG_list){
    
    print(date())
    
    # print cell type
    message(paste("Start working on ", celltype," in ",dir, " direction.",sep = ""))
    
    # subset object
    #obj <- subset(object, subset = cluster_celltype == celltype)
    obj <- object
    print(obj)
    print(table(obj$cluster_celltype))

    # gene expression normalization
    message("Noramlization in case not done before.")
    DefaultAssay(obj) <- "PC"
    obj <- SCTransform(obj,assay = "PC")

    #subset degs
    DEGs <- unique(DEG_list[DEG_list$dir == dir & DEG_list$celltype == celltype,]$gene)
    message(length(DEGs))

    # atac set normalization 
    DefaultAssay(obj) <- "CTpeaks"
    obj <- RunTFIDF(obj) %>% 
                FindTopFeatures(min.cutoff='q0') %>% 
                RunSVD()
    
    # first compute the GC content for each peak
    obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)

    # link peaks to genes
    message("Start calculating links.")
    
    obj <- LinkPeaks(
        object = obj,
        peak.assay = "CTpeaks",
        peak.slot = "data",
        expression.assay = "SCT",
        expression.slot = "data",
        genes.use = DEGs,
        method = "spearman",#c("pearson","spearman") # spearman may preferred.
        score_cutoff = 0.05, # set score to 0 to keep all 
        pvalue_cutoff = 0.05 # set p val = 1 to include all independent tests
    )

    # clean results and calculate adjust pvalues
    ## preform multiple testing correction
    message("Final processing the results.")
    res <- as.data.frame(Links(obj))
    res$p.adj <- p.adjust(res$pvalue,method = "BH")
    res <- res[res$score > 0.05 & res$p.adj < 0.05,]

    ## add annotataion information
    res$celltype <- celltype
    res$deg_dir <- dir
    res$comb <- paste(res$peak,res$gene,res$celltype,res$deg_dir, sep = "_")
    message(paste("Identified",length(res$comb),"links."))

    # reture
    return(res)
}

## actuall link analysis
ct <- c("Astrocyte","Excitatory","Inhibitory","Microglia","Oligodendrocyte","OPC")
direction <- c("pos","neg")

## running
final <- data.frame()
for (i in ct){
    for(j in direction){
        celltype <- i
        print(celltype)
        dir = j
        print(dir)
        
        res1 <- get_link(object=object, celltype=celltype, dir=dir, DEG_list=DEG_list)
        print(length(res1$gene))
        final <- rbind(final,res1)
    }
}

## CT peaks annotation files by MACS2
ctpeaks <- read.csv("./Results/CTpeaks_annotated.csv",row.names = 1)
colnames(ctpeaks) <- c("chr","start","end","width","strand","peak_called_in")
ctpeaks <- GRanges(ctpeaks)
ctpeaks

# pipeline to keep peaks only expressed in assigned celltype.  
get_clean_link <- function(links_df,ctpeaks){
    linked_peaks <- links_df

    # reformat the peak regions 
    linked_peaks$peak.start <-  str_split_fixed(linked_peaks$peak,"-",3)[,2]
    linked_peaks$peak.end <- str_split_fixed(linked_peaks$peak,"-",3)[,3]
    linked_peaks <- unique(linked_peaks)
    table(linked_peaks$celltype)

    final <- data.frame()

    ct <- c("Astrocyte","Excitatory","Inhibitory","Microglia","Oligodendrocyte","OPC")
    
    for(i in ct){
        # assign cell type
        celltype = i

        ## ct peaks that are not cell type AD-DEG linked peaks
        anno_temp <- ctpeaks[grep(celltype,ctpeaks$peak_called_in)]

        # get signals within specific cell type
        atac_peaks <- linked_peaks[linked_peaks$celltype == celltype,]

        link_temp <- GRanges(atac_peaks[,c("seqnames","peak.start","peak.end")])

        print(table(link_temp %in% anno_temp))
        # OUTPUT filtered results
        out <- atac_peaks[which(link_temp %in% anno_temp),]

        final <- rbind(final, out)
    }

    return(final)
}

# get peaks only shown in assigned celltype
#table(final$celltype, final$deg_dir)
links_clean <- get_clean_link(links_df = final, ctpeaks = ctpeaks)

## check information
table(links_clean$celltype,links_clean$deg_dir)

write.csv(links_clean, file = "./Results/LINK/linkpeaks_all_EC_1.23.csv",row.names = F)

################################################################################
# histogram plotting for linkage between peaks and genes in each brain region###
################################################################################
df <- table(paste(links_clean$gene,links_clean$celltype,sep = "_"))
#hist(df)
df <- as.data.frame(df)

avg.peaks <- round(median(df[,2]),digits = 0)
message(paste("Mean of number of peaks linked to DEGs: ",avg.peaks,".",sep = ""))

df$Freq <- ifelse(df$Freq >39,40,df$Freq)
# Represent it
p1 <- df %>%
  ggplot( aes(x=Freq)) +
    geom_histogram(fill="#825ca6ff", alpha=0.7, position = 'dodge',binwidth = 2,colour='#825ca6ff',size=1) +
    scale_fill_viridis(discrete=TRUE)+
    xlab("Linked peaks per gene")+
    ylab("Number of genes")+
    labs(fill="EC")+theme_classic()+ggtitle(paste("EC: ",avg.peaks,"linked peaks per DEG."))

pdf(file = "./Figures/LINK/Hist_link.peaks_EC_1.23.pdf",width = 5,height = 5)
p1
dev.off()

#EC "#825ca6ff"
#HIP "#3f78c199"
#PFC "#c25757ff"

## get annotation of the peaks linked to EOAD-DEG ##
## reference peaks bed Li et al. 2023 A comparative atlas of single-cell chrommatin accessibitlity in the human brain 
brain.peaks <- read.delim("./Data/cCREs.bed", header = F)
colnames(brain.peaks) <- c("chr","start","end","annotation")
brain.peaks <- GRanges(brain.peaks)
head(brain.peaks)

## ENCODE cCRE peaks #download from https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=encodeCcreCombined
encode_peaks <- import.bb("./Data/encodeCcreCombined.bb")
head(encode_peaks)

## add brain single cell eqtl information to the linked peaks
## get eqtl information
brain_eqlt <- read.csv("./Data/brain_sc_eqtl_sig.csv")
brain_eqlt <-  brain_eqlt[!is.na(brain_eqlt$chr),]
brain_eqlt$chr <- paste("chr",brain_eqlt$chr,sep = "")
dim(brain_eqlt)
#head(brain_eqlt)

# from expanded snps list, generate a GRanges format file
#snp_coords <- as.data.frame(sig[,c("chr","variant_pos","variant_pos","rs_id_dbSNP151_GRCh38p7","gene_name")])
snp_coords <- as.data.frame(brain_eqlt[,c("chr","pos","pos","SNP","symbol","cell_type","dist_TSS","beta")])

colnames(snp_coords) <- c("chr","start","end","RSID","gene_name","cell_type","dist_TSS","beta")
snp_coords <- GRanges(snp_coords)
head(snp_coords)


get_peak_annotation <- function(linked_peaks = linked_peaks, encode_peaks = encode_peaks, brain.peaks = brain.peaks){
    linked_peaks = linked_peaks
    ## change the format of linked peaks to GRanges
    #linked_peaks$peak.start <-  str_split_fixed(linked_peaks$peak,"-",3)[,2]
    #linked_peaks$peak.end <- str_split_fixed(linked_peaks$peak,"-",3)[,3]
    linked_gr <- linked_peaks[,c(1,14,15,4:13)]
    #linked_gr <- linked_peaks[,c(1,14,15,4:13,16)]
    df_granges <- GRanges(linked_gr)

    ## find overlap with encode and brain cCRE
    overlap_encode <- findOverlaps(df_granges,encode_peaks)
    overlap_brain_ccre <- findOverlaps(df_granges,brain.peaks)
    overlap_eqtl <- findOverlaps(df_granges,snp_coords)

    ## show results
    print("EOAD-DEG linked peaks overlap with ENCODE cis-CRE")
    print(length(unique(queryHits(overlap_encode))))
    print(length(unique(queryHits(overlap_encode)))/length(df_granges))

    print("EOAD-DEG linked peaks overlap with brain cis-CRE atlas")
    print(length(unique(queryHits(overlap_brain_ccre))))
    print(length(unique(queryHits(overlap_brain_ccre)))/length(df_granges))

    print("EOAD-DEG linked peaks overlap with single cell brain eqtl")
    # ATAC-seq peaks contain SNPs
    peaks_contain_snps <- linked_peaks[queryHits(overlap_eqtl),]
    peaks_contain_snps$RSID <- snp_coords[subjectHits(overlap_eqtl)]$RSID
    peaks_contain_snps$eGene <- snp_coords[subjectHits(overlap_eqtl)]$gene_name
    peaks_contain_snps$cell_type <- snp_coords[subjectHits(overlap_eqtl)]$cell_type
    peaks_contain_snps$dist_TSS <- snp_coords[subjectHits(overlap_eqtl)]$dist_TSS
    peaks_contain_snps$beta <- snp_coords[subjectHits(overlap_eqtl)]$beta
    peaks_contain_snps <- peaks_contain_snps[peaks_contain_snps$gene == peaks_contain_snps$eGene & peaks_contain_snps$celltype == peaks_contain_snps$cell_type,]

    print(length(unique(rownames(peaks_contain_snps))))
    print(length(unique(rownames(peaks_contain_snps)))/length(df_granges))

    ## add annotation if in ENCODE or brain cCRE atlas
    linked_peaks$in_encode <- FALSE
    linked_peaks[unique(queryHits(overlap_encode)),]$in_encode <- TRUE
    linked_peaks$in_brain_cCRE <- FALSE
    linked_peaks[unique(queryHits(overlap_brain_ccre)),]$in_brain_cCRE <- TRUE
    linked_peaks$in_brain_sc_eqlt <- ifelse(linked_peaks$comb %in% peaks_contain_snps$comb, T, F)

    ## get annotation from ENCODE
    q_df <- linked_peaks[queryHits(overlap_encode),]
    s_df <- as.data.frame(encode_peaks[subjectHits(overlap_encode)])

    # combine information 
    c_df <- cbind(q_df,s_df)
    c_df$encodeLabel <- factor(c_df$encodeLabel,levels = c("PLS","pELS","DNase-H3K4me3","dELS","CTCF-only","Other"))
    c_df <- c_df[order(c_df$encodeLabel,c_df$celltype,c_df$deg_dir,c_df$seqnames,c_df$start),]

    # based on the priority, keep the first one
    f_df <- c_df[!duplicated(c_df$comb),]

    ### add annotation information to original data
    id <- match(linked_peaks$comb,f_df$comb)

    linked_peaks$encodeLabel <- f_df[id,]$encodeLabel
    linked_peaks$ucscLabel <- f_df[id,]$ucscLabel
    linked_peaks$description <- f_df[id,]$description
    ## change na to other
    linked_peaks[is.na(linked_peaks$encodeLabel),]$encodeLabel <- "Other"
    linked_peaks[is.na(linked_peaks$ucscLabel),]$ucscLabel <- "Other"

    ## adding RSID information for brain eqtl
    # add rsid information for egene
    linked_peaks$eqtl_RSID <- NA
    linked_peaks$dist_TSS <- NA
    linked_peaks$beta <- NA
    id <- match(peaks_contain_snps$comb,linked_peaks$comb)
    linked_peaks[id,]$eqtl_RSID <- peaks_contain_snps$RSID
    linked_peaks[id,]$dist_TSS <- peaks_contain_snps$dist_TSS
    linked_peaks[id,]$beta <- peaks_contain_snps$beta

    return(linked_peaks)
}

### check if disease specific peaks are those novel peaks
linked_peaks <- read.csv("./Results/LINK/linkpeaks_all_EC_1.23.csv",header = 1)
linked_peaks <- linked_peaks[,-11]

p.adj <- p.adjust(linked_peaks$pvalue,method = "BH")
table(linked_peaks$p.adj < 0.05)
linked_peaks_annotated <- get_peak_annotation(linked_peaks = linked_peaks,encode_peaks = encode_peaks, brain.peaks = brain.peaks)
#pfc_linked_peaks_annotated[pfc_linked_peaks_annotated$in_brain_sc_eqlt,]
write.csv(linked_peaks_annotated, file = "./Results/LINK/EC_linkpeaks_all_annotated_1.23.csv", row.names = F)