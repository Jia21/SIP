#----------------------------------------------------------------------------------------------#
# Figure S5: gene expression of shared SIPs
# creates 
# requires shared_genes_allInfo.txt, rna_blood_cell.tsv, GTEx_mat_tissue_median
#----------------------------------------------------------------------------------------------#
library(data.table)
library(tidyr)
library(ggplot2)
library(ggpubr)

setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")
# setwd("/proj/yunligrp/users/lagler/SIP/data/")

shared <- fread("Shared_genes_allInfo.txt", drop=c(2,4:6))
uniqueN(shared,  "ENSEMBL_GENEID") #192 with ensembl id
uniqueN(shared,  "baitName") #234 with gene name
sharedlong <- melt(shared, id.vars = 1:2, 
                   variable.name = "tissue", value.name="TPM")
colnames(sharedlong)[[2]] <- "geneID"; colnames(sharedlong)[[1]] <- "geneName"
# only keep genes with expression info
sharedlong <- sharedlong[is.na(TPM)==F]

# hpa <- fread("rna_blood_cell.tsv", select=1:4)
# colnames(hpa) <- c("geneID", "geneName", "tissue", "TPM")
# #hpa2 <- hpa[geneName %in% shared$baitName | geneID %in% shared$ENSEMBL_GENEID]
# hpa2 <- hpa[geneName %in% shared$baitName &
#               tissue %in% c("classical monocyte", "intermediate monocyte",
#                             "naive CD4 T-cell", "neutrophil")]
# uniqueN(hpa2, "geneID"); uniqueN(hpa2, "geneName") #162

# /proj/yunligrp/users/jwen/SIP/post/GTEx_mat_tissue_median
gtex <- fread("GTEx_mat_tissue_median", drop=12)
colnames(gtex)[[1]] <- "geneID"; colnames(gtex)[[2]] <- "geneName"
#gtex2 <- gtex[geneName %in% shared$baitName | geneID %in% shared$ENSEMBL_GENEID]
gtex2 <- gtex[geneName %in% unique(sharedlong$geneName)]
uniqueN(gtex2, "geneName") #184
uniqueN(gtex2, "geneID") #185 
# one duplicate (SNORA79), remove the row of zeros
gtex2 <- gtex2[!(geneName=="SNORA79" & Thyroid==0)]

gtexlong <- melt(gtex2, id.vars = 1:2, 
                 variable.name = "tissue", value.name="TPM")

# combine data
exprdata <- rbind(gtexlong, sharedlong)
# unnecessary to do median since unique, but don't want to modify now. doesn't change anything
exprdata2 <- exprdata[, .("medTPM"=median(TPM, na.rm=T)), by=c("geneName", "tissue")]

setDT(exprdata2)[, color := ifelse(tissue=="Whole Blood", "Whole Blood",
                                  ifelse(tissue %in% c("Ery", "MacMon", "MK",
                                                       "nCD4", "Neu"), "Blood Cell", "Other Tissue"))]

# average by gene
avgGeneOther <- exprdata2[color!="Blood Cell",][, .("meanExpr"=mean(medTPM)), by="geneName"]
setDT(avgGeneOther)[, TopPerOther := cut(meanExpr, 
                               quantile(meanExpr, seq(.1, 1, .1)), 
                               labels = seq(90,10,-10),
                               include.lowest = F)]
avgGeneBlood <- exprdata2[color=="Blood Cell",][, .("meanExpr"=mean(medTPM)), by="geneName"]
setDT(avgGeneBlood)[, TopPerBlood := cut(meanExpr, 
                                         quantile(meanExpr, seq(.1, 1, .1), na.rm=T), 
                                         labels = seq(90,10,-10),
                                         include.lowest = F)]
#avgGeneOther[is.na(TopPercent)]$TopPercent <- 100

# add back to expression
exprdata3 <- merge(exprdata2, avgGeneOther[,c(1,3)], by="geneName")
exprdata3 <- merge(exprdata3, avgGeneBlood[,c(1,3)], by="geneName")
exprdata3$TopPerOther <- factor(exprdata3$TopPerOther, levels=seq(10,90,10))
exprdata3$TopPerBlood <- factor(exprdata3$TopPerBlood, levels=seq(10,90,10))

#----------------------------------------------------------------------------------------------#
# by GTEx expression
ggplot(exprdata3, aes(x=tissue, y=log10(medTPM+1), color=color)) +
  geom_violin()+
  geom_point(alpha=.6, size=1) +
  facet_wrap(~TopPerOther, ncol=1, scales="free") + 
  scale_color_manual(values=c("red", "black", "blue")) +  
  xlab("Tissue/Cell Type") + ylab("log10(TPM + 1)") +
  theme_light() + 
  theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1),
        legend.title = element_blank(), legend.position = "top", 
        plot.margin=unit(c(5.5,5.5,5.5,22),"points"))
ggsave("../figures/sharedSIPs_meanOther.pdf",
       width=15, height=70, unit="in", limitsize = F)

# by blood cell expression
ggplot(exprdata3, aes(x=tissue, y=log10(medTPM+1), color=color)) +
  geom_violin()+
  geom_point(alpha=.6, size=1) +
  facet_wrap(~TopPerBlood, ncol=1, scales="free") + 
  scale_color_manual(values=c("red", "black", "blue")) +  
  xlab("Tissue/Cell Type") + ylab("log10(TPM + 1)") +
  theme_light() + 
  theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1),
        legend.title = element_blank(), legend.position = "top", 
        plot.margin=unit(c(5.5,5.5,5.5,22),"points"))
ggsave("../figures/sharedSIPs_meanBlood.pdf",
       width=15, height=70, unit="in", limitsize = F)

#----------------------------------------------------------------------------------------------#
# Figures for Supplemental: only top 10 and 20

mkplot <- function(plotdata){
  ggplot(plotdata, aes(x=tissue, y=log10(medTPM+1), color=color)) +
    geom_violin()+
    geom_point(alpha=.6, size=.75) +
    scale_color_manual(values=c("red", "black", "blue")) +  
    xlab("Tissue/Cell Type") + ylab("log10(TPM + 1)") +
    theme_light() + 
    theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1),
          legend.title = element_blank(), legend.position = "top", 
          plot.margin=unit(c(5.5,5.5,5.5,22),"points"))
}  

ggarrange(
  mkplot(exprdata3[TopPerBlood==10])+ylim(0,4)+ggtitle("Top 10% of Blood Cells"),
  mkplot(exprdata3[TopPerBlood==20])+ylim(0,4)+ggtitle("Top 20% of Blood Cells"),
  nrow=2, labels=letters[1:2], align="hv", common.legend = T,
  legend = "top"
)
ggsave("../figures/sharedSIPs_topBlood.png", width=10, height=11, unit="in")

ggarrange(
  mkplot(exprdata3[TopPerOther==10])+ylim(0,4)+ggtitle("Top 10% of Other Tissues"),
  mkplot(exprdata3[TopPerOther==20])+ylim(0,4)+ggtitle("Top 20% of Other Tissues"),
  nrow=2, labels=letters[3:4], align="hv", common.legend = T,
  legend = "top"
)
ggsave("../figures/sharedSIPs_topOther.png", width=10, height=11, unit="in")


#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
genes <- unique(exprdata2$geneName)
ggplot(exprdata2, aes(x=tissue, y=log10(medTPM+1), color=color)) +
  geom_violin()+
  geom_point(alpha=.3) + 
  scale_color_manual(values=c("red", "black", "blue")) +  
  xlab("Tissue/Cell Type") + ylab("log10(TPM + 1)") +
  theme_light() + 
  theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1),
        legend.title = element_blank(), legend.position = "top", 
        plot.margin=unit(c(5.5,5.5,5.5,22),"points"))
ggsave("../figures/sharedSIPs.pdf", width=20, height=8, unit="in")


exprdata[TPM > 7500]
exprdata[geneName=="RPS27"]

max(exprdata2[tissue=="Whole Blood", medTPM])
exprdata2[medTPM==788] # RNA expression in bone marrow & lymphoid tissues & blood
hist(log10(exprdata2[geneName=="RHOG"]$medTPM+1))
