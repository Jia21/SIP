#----------------------------------------------------------------------------------------------#
# Figure 1: CDF plots of cell type specific genes versus shared genes
# requires CellTypeSpecificGene_*, shared.txt, union_ensgID.txt, and fraser.rda
#----------------------------------------------------------------------------------------------#
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr) #filter
setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")
#setwd("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/")
#----------------------------------------------------------------------------------------------#
# read in pcHiC data and shared gene info
#----------------------------------------------------------------------------------------------#
# pcHiC data
# formatted for intra-chromosomal interactions and distance < 2MB
load("fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]

# shared genes (entropy < 0.1)
shared <- fread("shared.txt", drop=1) #file generated with extra first column (Jia)
uniqueN(shared, by="ENSEMBL_GENEID") # 1559 (may include chrs outside of 1-22)

# add chr, start, end positions to merge with pcHiC data
union.id <- fread("union_ensgID.txt", drop=5)
union.id <- union.id[chr %in% 1:22]
shared2 <- merge(shared, union.id, by.x="ENSEMBL_GENEID", by.y="id")
uniqueN(shared2, by="ENSEMBL_GENEID") # 1476
setDT(shared2)[, type := "shared"]

#----------------------------------------------------------------------------------------------#
# Function returns empirical CDF plot and number of shared/specific genes 
#----------------------------------------------------------------------------------------------#
mkPlot <- function(cell, cellname){
  # specific genes by CT (entropy > 0.5, expression > 1)
  specific <- fread(paste0("CellTypeSpecificGene_", cell, ".txt"))
  #uniqueN(specific, by="ENSEMBL_GENEID") # 821
  
  # add chr, start, end positions to merge with pcHiC data
  specific2 <- merge(specific, union.id, by.x="ENSEMBL_GENEID", by.y="id")
  #uniqueN(specific2, by="ENSEMBL_GENEID") # 777
  setDT(specific2)[, type := "specific"]
  
  # combine shared and specific genes
  allgenes <- rbind(shared2[, c("ENSEMBL_GENEID", "BaitID", "type")],
                    specific2[, c("ENSEMBL_GENEID", "BaitID", "type")])
  allgenes <- unique(allgenes)
  #uniqueN(allgenes, by="ENSEMBL_GENEID") # 2253 (specific + shared)
  
  
  # filter for significant interactions based on CT
  hic <- fraser %>% filter(!!as.name(cell) >=5) %>% data.table
  # merge with pcHiC data to count interactions
  hic2 <- merge(hic, allgenes, by.x="baitID", by.y="BaitID") 
  #uniqueN(hic2, by="ENSEMBL_GENEID") # 2003 a gene could correspond to 0 significant interactions
  
  nShared <- uniqueN(hic2[type=="shared"], by="ENSEMBL_GENEID") # 1322
  nSpecific <- uniqueN(hic2[type=="specific"], by="ENSEMBL_GENEID") # 681
  
  
  # number of interactions per gene
  nInter <- hic2[, .N, by="ENSEMBL_GENEID"]
  
  # add back to allgenes data
  allgenesN <- merge(allgenes, nInter)
  setorder(allgenesN, -type, N)
  
  # t test for difference in # of interactions
  t <- t.test(allgenesN[type=="specific", N], allgenesN[type=="shared", N])
  est.sp <- format(as.numeric(t$estimate[[1]]), digits = 3)
  est.sh <- format(as.numeric(t$estimate[[2]]), digits = 3)
  pval <- format(t$p.value, digits=2)
  
  
  # empirical CDF
  ecdf <-  ggplot(allgenesN, aes(N, color=type)) +
    stat_ecdf(pad=F, size=.7) + 
    scale_color_manual(values=c("blue", "red"), labels=c("Shared Genes", "Specific Genes")) +
    xlab("# of Interactions") + ylab("CDF") + ggtitle(cellname) + 
    annotate("label", x=max(allgenesN$N)*.6, y=.5, size=3.5,
             label=paste0("Avg. Specific = ", est.sp,
                          "\nAvg. Shared = ", est.sh,
                          "\np = ", pval))+
    theme_bw() + theme(legend.title=element_blank(), legend.text=element_text(size=12))
 
  return(list("plot"=ecdf,
              "nGenes"= c(cellname, format(nSpecific, big.mark=","), format(nShared, big.mark=","))))
  
}

#----------------------------------------------------------------------------------------------#
# call function for each cell type group
#----------------------------------------------------------------------------------------------#
# make table with number of specific/shared genes per cell type
nGenes <- rbind(mkPlot("Ery", "Erythrocyte")[[2]],
                mkPlot("Mac_Mon", "Macrophage/Monocyte")[[2]],
                mkPlot("MK", "Megakaryocyte")[[2]],
                mkPlot("nCD4", "Naive CD4 T-cell")[[2]],
                mkPlot("Neu", "Neutrophil")[[2]])
colnames(nGenes) <- c("Cell Type", "# Specific", "# Shared")
ggtable <- ggtexttable(nGenes[,1:2], theme=ttheme("lBlack", base_size=10))

# arrange all eCDF plots and table
ggarrange(mkPlot("Ery", "Erythrocyte")[[1]],
          mkPlot("Mac_Mon", "Macrophage/Monocyte")[[1]],
          mkPlot("MK", "Megakaryocyte")[[1]],
          mkPlot("nCD4", "Naive CD4 T-cell")[[1]],
          mkPlot("Neu", "Neutrophil")[[1]], 
          ggtable,
          nrow=2, ncol=3, labels=letters[1:6],
          common.legend = T)
ggsave("../figures/FigureS1_eCDF.pdf", width=10, height=7)
ggsave("../figures/FigureS1_eCDF.png", width=10, height=7)
#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#

