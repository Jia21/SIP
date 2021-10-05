##########################################################################
### Revised from /proj/yunligrp/users/lagler/SIP/R/hugin2.R ###
###                   By Jia Wen 09/2021                               ###
##########################################################################

setwd("/proj/yunligrp/users/jwen/SIP/revise/1_4")

#----------------------------------------------------------------------------------------------#
# Figure 4: SIP subnetworks 
# creates figures similar to Hugin2 with SNPs, genes, and pcHiC interactions
# requires fraser.rda, atac seq data, and SIPsubnetworks.xlsx
#----------------------------------------------------------------------------------------------#
library(data.table)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(xlsx)
library(viridis)
library(scales)
library(ggrepel)

load("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[[31]] <- "MacMon"

expr <- fread("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/GeneExpressionMatrix_combine_Mac_Mon.txt") 

mk <- read.xlsx2("/proj/yunligrp/users/lagler/SIP/data/SIPsubnetworks.xlsx", sheetName = "MK_PLT")%>%
  subset(., baitID=="492452") %>% merge(., fraser[,c("baitID", "oeID", "MK", "dist")],
            by=c("baitID", "oeID"), all.x=T, all.y=F) %>% data.table()
mk$oeStart <- as.numeric(as.character(mk$oeStart))
mk$oeEnd <- as.numeric(as.character(mk$oeEnd))
mk$baitStart <- as.numeric(as.character(mk$baitStart))
mk$baitEnd <- as.numeric(as.character(mk$baitEnd))
mk$MK <- as.numeric(as.character(mk$MK))
mk$POS_.hg19. <- as.numeric(as.character(mk$POS_.hg19.))

# interaction level data
mkU <- unique(mk, by="oeID")
setDT(mkU)[, xstart := (baitStart + baitEnd)/2]
setDT(mkU)[, xend := (oeStart + oeEnd)/2]
setorder(mkU, dist)
mkU$curvature <- c(-.15, -.1, -.05)

gene <- data.table("name"=c("CLCN2", "POLR2H", "CHRD", "THPO", "EPHB3"),
                   "start"=c(184079439, 184079506, 184097861, 184095932, 184279572),
                   "end"=c(184063973, 184086384, 184108524, 184089723, 184300197),
                   "ensemblID"=c("ENSG00000114859", "ENSG00000163882",
                                 "ENSG00000090539", "ENSG00000090534", "ENSG00000182580"),
                   "expr"=c(0.6218563,25.56251,0.0218635,0.001710853, 0.0005911997))

# atac-seq data
# /proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/Mega_peaks.narrowPeak.gz
atacmk <- fread("/proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/Mega_peaks.narrowPeak.gz", select = c(1:3))
colnames(atacmk) <- c("chr", "start", "end")
atacmk <- atacmk[chr == "chr3"]
setkey(mkU, oeStart, oeEnd)
atacmk_ovlp <- foverlaps(atacmk, mkU, by.x=c("start", "end"), nomatch=NULL) %>%
  unique(., by=c("start", "end"))


chicago <- 
  ggplot() +
  # shade PIRs
  geom_rect(data=mkU,
            aes(xmin=oeStart, xmax=oeEnd, ymin=4, ymax=max(mkU$MK)), alpha=.1) +
  # shade SIP bait
  geom_rect(data=unique(mk, by="baitID"),
            aes(xmin=baitStart, xmax=baitEnd, ymin=4, ymax=max(mkU$MK)), alpha=.4) +
  # significance line at y=5
  geom_hline(yintercept=5, size=.3, linetype="dashed", color="blue", alpha=.4) + 
  # draw chicago scores over interactions
  geom_segment(data=mk, aes(x=oeStart, xend=oeStart, y=4, yend=MK), color="blue") + 
  geom_segment(data=mk, aes(x=oeStart, xend=oeEnd, y=MK, yend=MK), color="blue") + 
  geom_segment(data=mk, aes(x=oeEnd, xend=oeEnd, y=4, yend=MK), color="blue") +
  # draw interaction lines
  lapply(split(mkU, 1:3), function(dat) {
    geom_curve(data=dat, aes(x=xstart, xend=xend, y=4, yend=4),
               curvature = dat$curvature, color="purple") }) +
  # add line at y=0 and modify theme
  geom_hline(yintercept=4, size=.25)+
  #xlim(min(mk$oeStart)-7500, max(mk$baitEnd)+20000) +
  scale_x_continuous(breaks=pretty_breaks(n=8), limits=c(min(mk$oeStart)-7500, max(mk$baitEnd)+20000)) +
  scale_y_continuous(breaks=seq(5,15,5), limits=c(-.5, max(mkU$MK))) +
  xlab("Position") + ylab("CHiCAGO Score") + 
  theme_bw() + theme(plot.margin = unit(c(-.5, 5.5, 15, 5.5), "points"),
                     axis.title.y = element_text(size=7),
                     axis.title.x = element_text(size=9),
                     axis.text.x = element_text(size=7))

# rsIDs
setorder(mk, POS_.hg19.)
rsid <- ggplot(data=mk) + 
  geom_segment(aes(x=POS_.hg19., xend=POS_.hg19., y=c(0:7), yend=c(1:8)))+
  geom_text_repel(aes(x=POS_.hg19., y=c(0.5:7.5), label=rsID),vjust = 1, nudge_y = -2,max.overlaps = Inf) +
  xlim(min(mk$oeStart)-7500, max(mk$baitEnd)+20000) +
  theme_void() + theme(panel.background=element_rect(colour="black", fill=NA),
                       axis.title.y = element_text(angle=90)) + 
  ylab("Relevant Variants")

  geneexpr <- 
  ggplot(data=gene) + 
  geom_segment(aes(x=start, xend=end, y=c(1,1.5,1.5,1,1), yend=c(1,1.5,1.5,1,1), color=expr),
               lineend="butt", linejoin="mitre", size=6,
               arrow=arrow(length=unit(.07, "inches"))) +
  scale_color_viridis(option="C") + labs(color="Expression")+
  guides(color = guide_colorbar(title.position = "top")) + 
  geom_text(aes(x=start, y=c(1,1.5,1.5,1,1)+.15, label=name), size=3) + 
  xlim(min(mk$oeStart)-7500, max(mk$baitEnd)+20000) + ylim(0.9,1.75)+
  theme_void() + theme(panel.background=element_rect(colour="black", fill=NA),
                       plot.margin = unit(c(-.04, 0, 0, 0), "lines"),
                       legend.position = c(.92,.75),
                       legend.direction = "horizontal",
                       axis.title.y = element_text(angle=90, size=7),
                       legend.title = element_text(size=6),
                       legend.key.size = unit(.2, 'in'),
                       legend.text = element_text(size=6)) +
  ylab("Genes")

# atac-seq
atac <- ggplot(data=atacmk_ovlp)+
  geom_segment(aes(x=start, xend=end, y=rep(1,7), yend=rep(1,7)), size=4) + 
  xlim(min(mk$oeStart)-7500, max(mk$baitEnd)+20000) +
  theme_void() + theme(panel.background=element_rect(colour="black", fill=NA),
                       plot.margin = unit(c(-.04, 0, 0, 0), "lines"),
                       axis.title.y = element_text(angle=90, size=7)) +
  ylab("ATAC")

# arrange all plots and save
mkplot <- ggarrange(rsid, geneexpr, atac, chicago,
          nrow=4, align="v", heights = c(.09,.225, .06, .625))
ggsave("Fig_4a_HuginMKPLT_rsID_1.png", height=6, width=12)


#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
# EXAMPLE 2: NCD4 LYM
#----------------------------------------------------------------------------------------------#
# sip subnetwork data
# setwd("/proj/yunligrp/users/lagler/SIP/data/")
nc <- read.xlsx2("/proj/yunligrp/users/lagler/SIP/data/SIPsubnetworks.xlsx", sheetName = "NCD4_LYM")%>%
  subset(., baitID=="139778") %>% merge(., fraser[,c("baitID", "oeID", "nCD4", "dist")],
                                        by=c("baitID", "oeID"), all.x=T, all.y=F) %>% data.table()
nc$oeStart <- as.numeric(as.character(nc$oeStart))
nc$oeEnd <- as.numeric(as.character(nc$oeEnd))
nc$baitStart <- as.numeric(as.character(nc$baitStart))
nc$baitEnd <- as.numeric(as.character(nc$baitEnd))
nc$nc <- as.numeric(as.character(nc$nc))
nc$POS_.hg19. <- as.numeric(as.character(nc$POS_.hg19.))

# interaction level data is the same here
setDT(nc)[, xstart := (baitStart + baitEnd)/2]
setDT(nc)[, xend := (oeStart + oeEnd)/2]
setorder(nc, dist)
nc$curvature <- c(-.125, -.1, -.1)

# no other end genes listed. only bait gene
gene2 <- data.table("name"=c("ETS1"),
                   "start"=c(128457453),
                   "end"=c(128328656),
                   "ensemblID"=c("ENSG00000134954"),
                   "expr"=c(325.4768))

# atac-seq data
# /proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/CD4_peaks.narrowPeak.gz
atacnc <- fread("/proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/CD4_peaks.narrowPeak.gz", select = c(1:3))
colnames(atacnc) <- c("chr", "start", "end")
atacnc <- atacnc[chr == "chr11"]
setkey(nc, oeStart, oeEnd)
atacnc_ovlp <- foverlaps(atacnc, nc, by.x=c("start", "end"), nomatch=NULL) %>%
  unique(., by=c("start", "end"))

#----------------------------------------------------------------------------------------------#

# chicago score interactions
max <- max(nc$nCD4)
chicago2 <- 
  ggplot() +
  # shade PIRs
  geom_rect(data=nc,
            aes(xmin=oeStart, xmax=oeEnd, ymin=4, ymax=max), alpha=.1) +
  # shade SIP bait
  geom_rect(data=unique(nc, by="baitID"),
            aes(xmin=baitStart, xmax=baitEnd, ymin=4, ymax=max), alpha=.4) +
  # significance line at y=5
  geom_hline(yintercept=5, size=.3, linetype="dashed", color="blue", alpha=.4) + 
  # draw chicago scores over interactions
  geom_segment(data=nc, aes(x=oeStart, xend=oeStart, y=4, yend=nCD4), color="blue") + 
  geom_segment(data=nc, aes(x=oeStart, xend=oeEnd, y=nCD4, yend=nCD4), color="blue") + 
  geom_segment(data=nc, aes(x=oeEnd, xend=oeEnd, y=4, yend=nCD4), color="blue") +
  # draw interaction lines
  lapply(split(nc[1:2,], 1:2), function(dat) {
    geom_curve(data=dat, aes(x=xstart, xend=xend, y=4, yend=4),
               curvature = dat$curvature, color="purple") }) +
    geom_curve(data=nc[3,], aes(x=xend, xend=xstart, y=4, yend=4),
             curvature = nc[3,curvature], color="purple") +
  #add line at y=0 and theme formatting
  geom_hline(yintercept=4, size=.25)+
  #xlim(min(nc$oeStart), max(nc$oeEnd)) + 
  xlab("Position") + ylab("CHiCAGO Score") + 
  scale_y_continuous(breaks=seq(4,6,1), limits = c(3.4,max)) + 
  scale_x_continuous(breaks=pretty_breaks(n=7), limits=c(min(nc$oeStart), max(nc$oeEnd))) +
  theme_bw() + theme(plot.margin = unit(c(-.5, 5.5, 5.5, 5.5), "points"),
                     axis.title.y = element_text(size=7),
                     axis.title.x = element_text(size=9),
                     axis.text.x = element_text(size=7))

# rsIDs
setorder(nc, POS_.hg19.)
rsid2 <- ggplot(data=nc) + 
  geom_segment(aes(x=POS_.hg19., xend=POS_.hg19., y=0, yend=.5))+
  geom_label(aes(x=POS_.hg19., y=c(0.5, 0.5, 0.5), label=rsID)) + 
  xlim(min(nc$oeStart), max(nc$oeEnd)) + 
  ylim(-.5,1.5)+
  theme_void() + theme(panel.background=element_rect(colour="black", fill=NA),
                       axis.title.y = element_text(angle=90, size=7),
                       plot.margin = unit(c(11, 5.5, 0, 5.5), "points")) + 
  ylab("SNPs")

# gene expression
geneexpr2 <- 
  ggplot(data=gene2) + 
  geom_segment(aes(x=start, xend=end, y=1, yend=1, color=expr),
               lineend="butt", linejoin="mitre", size=6,
               arrow=arrow(length=unit(.07, "inches"))) +
  scale_color_viridis(option="C", limits=c(0,350)) + 
  labs(color="Expression")+ guides(color = guide_colorbar(title.position = "top")) + 
  geom_text(aes(x=start, y=1.055, label=name), size=3.25) +
  xlim(min(nc$oeStart), max(nc$oeEnd)) + ylim(0.95,1.15) + 
  theme_void() + theme(panel.background=element_rect(colour="black", fill=NA),
                       plot.margin = unit(c(-.04, 0, 0, 0), "lines"),
                       legend.position = c(.92,0.65),
                       legend.direction = "horizontal",
                       axis.title.y = element_text(angle=90, size=7),
                       legend.title = element_text(size=6),
                       legend.key.size = unit(.2, 'in'),
                       legend.text = element_text(size=6)) +
  ylab("Genes")

# atac-seq
atac2 <- ggplot(data=atacnc_ovlp)+
  geom_segment(aes(x=start, xend=end, y=rep(1,11), yend=rep(1,11)), size=4) + 
  xlim(min(nc$oeStart), max(nc$oeEnd)) +
  theme_void() + theme(panel.background=element_rect(colour="black", fill=NA),
                       plot.margin = unit(c(-.04, 0, 0, 0), "lines"),
                       axis.title.y = element_text(angle=90, size=7)) +
  ylab("ATAC")

# arrange all plots and save
ncd4Plot <- ggarrange(rsid2, geneexpr2, atac2, chicago2,
          nrow=4, align="v", heights = c(.12,.2,.08,.6))
ggsave("Fig4b_HuginNCD4LYM_rsID.png", height=5, width=12)

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
# both plots
ggarrange(mkplot, ncd4Plot, nrow=2,
        labels=letters[1:2], vjust=1, hjust=0, 
        heights = c(.525, .475), align="hv")
ggsave("Figure4_subnetworks2_rsID.png", height=8, width=8)
ggsave("Figure4_subnetworks2_rsID.pdf", height=8, width=8)

### overlap with ATAC-seq ###
module add bedtools
cell=MK
zcat /proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/Mega_peaks.narrowPeak.gz |cut -f1-3 | sed 's/^chr//g' > ${cell}_ATAC.bed
bedtools intersect -wao -a ${cell}_GWASvars.bed -b ${cell}_ATAC.bed > ${cell}_ATAC_olap

cell=Ery
zcat /proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/${cell}_peaks.narrowPeak.gz |cut -f1-3 | sed 's/^chr//g' > ${cell}_ATAC.bed
bedtools intersect -wao -a ${cell}_GWASvars.bed -b ${cell}_ATAC.bed > ${cell}_ATAC_olap

cell=CD4
zcat /proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/${cell}_peaks.narrowPeak.gz |cut -f1-3 | sed 's/^chr//g' > ${cell}_ATAC.bed
bedtools intersect -wao -a ${cell}_GWASvars.bed -b ${cell}_ATAC.bed > ${cell}_ATAC_olap

cell=MacMon
zcat /proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/Mono_peaks.narrowPeak.gz |cut -f1-3 | sed 's/^chr//g' > ${cell}_ATAC.bed
bedtools intersect -wao -a ${cell}_GWASvars.bed -b ${cell}_ATAC.bed > ${cell}_ATAC_olap

# For Neu peaks, we used /proj/yunligrp/users/lagler/SIP/old/ATACseq/Neu_controlUnion_peaks.narrowPeak, but used in the ms finally #

#### For motifs overlaping ####
zcat /proj/yunligrp/users/quansun/annotation/vampire/motif_anal/ukbb_mbreaker_finemap_merged.txt.gz |cut -f3 |sort | uniq -c 
#   7573 baso
#    4547 baso_p
#   22192 eo
#   18184 eo_p
#   20489 hct
#   23364 hgb
#   16304 hlr
#   15056 hlr_p
#   10151 irf
#   22698 lymph
#   12490 lymph_p
#   22372 mch
#    8762 mchc
#   24456 mcv
#   23655 mono
#   12316 mono_p
#   21503 mpv
#   22601 mrv
#   20792 mscv
#   19523 neut
#   11390 neut_p
#   19118 pct
#   10131 pdw
#   28620 plt
#   27209 rbc
#   20823 rdw_cv
#   19166 ret
#   15567 ret_p
#   26154 wbc

Ery -> rdw_cv, hct, mch, hgb, mchc, rbc
MacMon -> mono, mono_p, wbc
MK -> plt, mpv
nCD4 <- "lymph", "lymph_p"
Neu <- "neut", "neut_p","wbc"

motifs=/proj/yunligrp/users/quansun/annotation/vampire/motif_anal/ukbb_mbreaker_finemap_merged.txt.gz
for cell in Ery MacMon MK nCD4 Neu; do 
    for phenotype in rdw_cv hct mch hgb mchc rbc mono mono_p wbc plt mpv lymph lymph_p neut neut_p; do
        echo ${cell}
        if [ ${cell} == "Ery" ];then
            if [ ${phenotype} == "rdw_cv" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi

            if [ ${phenotype} == "hct" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi

            if [ ${phenotype} == "mch" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi

            if [ ${phenotype} == "hgb" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi

            if [ ${phenotype} == "mchc" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi

            if [ ${phenotype} == "rbc" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi
        fi

        if [ ${cell} == "MacMon" ];then
            if [ ${phenotype} == "mono" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi

            if [ ${phenotype} == "mono_p" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi
            if [ ${phenotype} == "wbc" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi           
        fi

        if [ ${cell} == "MK" ];then
            if [ ${phenotype} == "plt" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi

            if [ ${phenotype} == "mpv" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi    
        fi

        if [ ${cell} == "nCD4" ];then
            if [ ${phenotype} == "lymph" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi

            if [ ${phenotype} == "lymph_p" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi    
        fi

        if [ ${cell} == "Neu" ];then
            if [ ${phenotype} == "neut" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi

            if [ ${phenotype} == "neut_p" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi    
            if [ ${phenotype} == "wbc" ]; then
                zcat ${motifs} | awk -v var=${phenotype} '$3==var' | cut -f1,13 >> ${cell}_motifs_vars
            fi 
        fi
    done
done

for cell in Ery MacMon MK nCD4 Neu; do
    sort ${cell}_motifs_vars |uniq > ${cell}_motifs_vars_uniq
done

### 
rm(list=ls())
library(tidyverse)
cell<-"Neu"
vars <- read.table(paste0(cell,"vars.bed"),header =T)
motifs<-read.table(paste0(cell,"_motifs_vars_uniq"),header =F)
vars_motifs <- left_join(vars,motifs,by=c("VariantID"="V1"))
vars$genename <- "NA"
for(i in 1:nrow(vars)){
    cat(i,"\t")
    tmp <- vars_motifs %>% filter(baitID== vars$baitID[i]) %>% filter(VariantID==vars$VariantID[i])
    vars$genename[i] <- paste(tmp$V2,collapse=";")
}
# vars_motifs_combine <- vars_motifs %>% group_by(baitID,VariantID) %>% summarise(paste(V2,collapse=";"))
fwrite(vars,paste0(cell,"_vars_motifs"),sep="\t")