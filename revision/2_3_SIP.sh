#!/bin/bash

cd /proj/yunligrp/users/jwen/SIP/revise

### Extract SIP and non-SIP regions
setwd("/proj/yunligrp/users/jwen/SIP/revise")
load("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[[31]] <- "MacMon"

## Load SIP results for 5 cell types ##
CT <- c("Ery","MK","MacMon","Neu","nCD4")
for (i in 1:length(CT)){
    cat(CT[i],"\n")
    sips <- read.table(paste0("/proj/yunligrp/users/lagler/SIP/data/SIP_",CT[i],".txt"),header = T) %>% filter(SIP==1)
    sip_baitID <- unique(sort(sips$baitID))

    fraser_ct <- fraser[eval(as.name(CT[i])) >= 5,]
    non_sip <- unique(setdiff(fraser_ct$baitID, sip_baitID))

    fraser_ct.bed <- fraser_ct %>% filter(baitID %in% sip_baitID) %>% select("baitChr","baitStart","baitEnd") %>% distinct()
    fraser_ct.bed$baitChr <- paste0("chr",fraser_ct.bed$baitChr)
    fraser_ct_nonsip.bed <- fraser_ct %>% filter(baitID %in% non_sip) %>% select("baitChr","baitStart","baitEnd") %>% distinct()
    fraser_ct_nonsip.bed$baitChr <- paste0("chr",fraser_ct_nonsip.bed$baitChr)
    write.table(fraser_ct.bed,paste0("./2_3/",CT[i],"_SIP.bed"),quote = F,sep="\t",col.name = F,row.names = F)
    write.table(fraser_ct_nonsip.bed,paste0("./2_3/",CT[i],"_nonSIP.bed"),quote = F,sep="\t",col.name = F,row.names = F)
}

### Extract the sequence ###
zcat /proj/yunligrp/users/quansun/annotation/vampire/motif_anal/ukbb_mbreaker_finemap_merged.txt.gz | \
     cut -f4-6,13,17 | sort | uniq | grep -w "GATA1" > /proj/yunligrp/users/jwen/SIP/revise/2_3/ukbb_mbreaker_finemap_merged_motifs_GATA1.bed

### generate the meme file format ###
wget http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme 
# GATA1 # 
sed -n '3013,3034P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > GATA1.meme
# KLF1 #
sed -n '5127,5143P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > KLF1.meme
# KLF6 #
sed -n '5200,5221P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > KLF6.meme
# MAFB #
sed -n '5482,5495P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > MAFB.meme
# ATF3 #
sed -n '408,421P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > ATF3.meme
# GATA2 #
sed -n '3051,3072P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > GATA2.meme
# FLI1 #
sed -n '2417,2437P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > FLI1.meme
# RUNX1 #
sed -n '9083,9096P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > RUNX1.meme
# LYL1 #
sed -n '5439,5455P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > LYL1.meme
# ERG #
sed -n '2062,2077P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > ERG.meme
# GFI1B #
sed -n '3266,3278P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > GFI1B.meme
# TAL1 #
sed -n '10289,10329P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > TAL1.meme
# EVI1 #
sed -n '2352,2370P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > EVI1.meme
# ETV6 #
sed -n '2328,2339P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > ETV6.meme
# GATA3 #
sed -n '3089,3102P' HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > GATA3.meme

module add bedtools

fastaFromBed -fi /proj/yunligrp/users/jwen/Reference/hg19.fa -bed /proj/yunligrp/users/jwen/SIP/revise/2_3/ukbb_mbreaker_finemap_merged_motifs_GATA1.bed -fo /proj/yunligrp/users/jwen/SIP/revise/2_3/ukbb_motifs_GATA1.fa
for cell in Ery MacMon Neu nCD4 MK;do
    echo ${cell}
    fastaFromBed -fi /proj/yunligrp/users/jwen/Reference/hg19.fa -bed ${cell}_SIP.bed -fo /proj/yunligrp/users/jwen/SIP/revise/2_3/${cell}_SIPs.fa
    fastaFromBed -fi /proj/yunligrp/users/jwen/Reference/hg19.fa -bed ${cell}_nonSIP.bed -fo /proj/yunligrp/users/jwen/SIP/revise/2_3/${cell}_nonSIPs.fa
done

module add meme/5.2.0
for cell in MacMon Ery Neu MK nCD4;do
    for m in GATA1 MAFB KLF1 KLF6 ATF3 GATA2 FLI1 RUNX1 ERG GFI1B TAL1 EVI1 ETV6 GATA3; do
        echo ${cell}, ${m}
        fimo -o /proj/yunligrp/users/jwen/SIP/revise/2_3/${cell}_SIP_${m} ${m}.meme ${cell}_SIPs.fa
        fimo -o /proj/yunligrp/users/jwen/SIP/revise/2_3/${cell}_nonSIP_${m} ${m}.meme ${cell}_nonSIPs.fa
    done
done

### Put together all the motif enrichment results ###
setwd("/proj/yunligrp/users/jwen/SIP/revise/2_3")
rm(list=ls())
cell <- c("Ery","MacMon","Neu","MK","nCD4")
motifs <- c("FLI1","KLF1","KLF6","GATA1", "MAFB", "ATF3", "GATA2", "RUNX1", "ERG", "GFI1B" ,"TAL1", "EVI1", "ETV6" ,"GATA3")

count_table <- NULL
for (c in 1:length(cell)){
    for(m in 1:length(motifs)){
        cat(cell[c],motifs[m],"\n")
        fa <- read.table(paste0("/proj/yunligrp/users/jwen/SIP/revise/2_3/",cell[c],"_SIPs.fa"),header = F)
        nosip_fa <- read.table(paste0("/proj/yunligrp/users/jwen/SIP/revise/2_3/",cell[c],"_nonSIPs.fa"),header = F)
        dat_sip <- read.table(paste0("/proj/yunligrp/users/jwen/SIP/revise/2_3/",cell[c],"_SIP_",motifs[m],"/fimo.tsv"),header = T,sep="\t",fill=T)
        dat_nonsip <- read.table(paste0("/proj/yunligrp/users/jwen/SIP/revise/2_3/",cell[c],"_nonSIP_",motifs[m],"/fimo.tsv"),header = T,sep="\t",fill=T)
        count1 <- length(unique(sort(dat_sip$sequence_name)))
        count2 <- length(unique(sort(dat_nonsip$sequence_name)))
        count_table <- rbind(count_table,c(cell[c],motifs[m],count1,nrow(fa)/2,count2,nrow(nosip_fa)/2))
    }
}
count_table <- data.frame(count_table)
count_table$X3 <- as.character(count_table$X3)
count_table$X4 <- as.character(count_table$X4)
count_table$X5 <- as.character(count_table$X5)
count_table$X6 <- as.character(count_table$X6)
count_table$p <- 0
count_table$or <- 0
count_table$lower <- 0
count_table$upper <- 0
for(i in 1:nrow(count_table)){
    test <- matrix(c(as.numeric(count_table[i,3]),(as.numeric(count_table[i,4])-as.numeric(count_table[i,3])),
                   as.numeric(count_table[i,5]),(as.numeric(count_table[i,6])-as.numeric(count_table[i,5]))),
                   nrow = 2, byrow=T,
                   dimnames = list(SIP = c("SIP", "nonSIP"),Motifs = c("enriched", "noenriched")))
    #  test1 <- matrix(c(as.numeric(count_table[i,3]),(as.numeric(count_table[i,4])),
    #                as.numeric(count_table[i,5]),(as.numeric(count_table[i,6]))),
    #                nrow = 2, byrow=T,
    #                dimnames = list(SIP = c("SIP", "nonSIP"),Motifs = c("enriched", "noenriched")))
    p <- fisher.test(test,alternative="greater")$p.value
    or <- fisher.test(test,alternative="greater")$estimate
    low <- fisher.test(test,alternative="greater")$conf.int[1]
    upp <- fisher.test(test,alternative="greater")$conf.int[2]
    count_table$p[i] <- p
    count_table$or[i] <- or
    count_table$lower[i] <- low
    count_table$upper[i] <- upp
}


library(tidyverse)
count_table_sub <- count_table %>% filter(X2 %in% c("FLI1","KLF1","KLF6","ERG","TAL1","EVI1","ETV6","GATA1","RUNX1","GATA2"))
count_table_sub$log <- -log10(as.numeric(as.character(count_table_sub$p)))/8
p <- ggplot(count_table_sub, aes(x = X2, y = or,ymin = or, ymax = or,label = formatC(p,format = "e", digits = 2))) +
     geom_pointrange(position = position_dodge(0.5),fatten=6) + 
     geom_hline(yintercept = 1, lty = 2) + 
     coord_flip() + ylim(0,11) +
     facet_wrap(.~X1) +
     xlab("TF") + ylab("Enrichment (95% CI)") +
     theme_bw() + 
     geom_text(vjust = -1.5,  hjust = 0.6, size = 4) +
     theme(axis.text=element_text(size=15),axis.title=element_text(size=15),strip.text.x = element_text(size = 15)) 
ggsave("motif_enrich.pdf", height=12, width=18)
ggsave("motif_enrich.png", height=12, width=18)

############################################################
#############
# 1. TADs #
#############
library(data.table)
library(tidyr) # separate
library(readxl) # read_xlsx
library(tidyverse)
tissues <- c("GM12878")
TADs <- NULL
for(i in 1:length(tissues)){
    cat(tissues[i],"\n")
    dat <- read_xlsx("/proj/yunligrp/users/jwen/TOPMed_SV/anno/mmc4.xlsx",col_names = F,
                   sheet=tissues[i]) %>% data.table()
    colnames(dat) <- c("chr","st","end")   
    TADs <- rbind(TADs,cbind(dat,tissues[i]))
}
TADs = format(TADs,scientific=F)
TADs <- gsub(" ","",TADs)
TADs <- data.frame(TADs)
write.table(TADs,"TADs_GM12878_hg19.bed",quote = F,sep="\t",col.names = F,row.names =F)

module add bedtools
sort -k1,1 -k2,2n TADs_GM12878_hg19.bed > TADs_GM12878_hg19_sort.bed
for cell in Ery MacMon Neu nCD4 MK; do
    echo ${cell}
    bedtools closest -a ${cell}_SIP.bed -b TADs_GM12878_hg19_sort.bed -D a > ${cell}_sip.dist.TADs
    bedtools closest -a ${cell}_nonSIP.bed -b TADs_GM12878_hg19_sort.bed -D a > ${cell}_nonsip.dist.TADs
done

setwd("/proj/yunligrp/users/jwen/SIP/revise/2_3")
library(ggplot2)
cell <- c("Ery","MacMon","Neu","MK","nCD4")
p <- NULL
for(i in 1:length(cell)){
    cat(cell[i],"\n")
    dat <- read.table(paste0(cell[i],"_sip.dist.TADs"),header = F)
    dat_nonsip <- read.table(paste0(cell[i],"_nonsip.dist.TADs"),header = F)
    dat$label <- "SIP"
    dat$kb <- dat$V8/1000
    dat_nonsip$label <- "non-SIP"
    dat_nonsip$kb <- dat_nonsip$V8/1000
    plotd <- rbind(dat,dat_nonsip)

    ggplot(plotd, aes(x=kb, color=label)) +
        geom_density() + scale_color_brewer(palette="Paired") + xlab("Distance from SIP/non-SIP to TADs (KB)") + ylab("Density") +
        scale_y_continuous(breaks = c(0, 5e-04, 1e-03),labels = c("0",expression(5^-4),expression(1^-3))) +  ggtitle(cell[i]) +
        theme_bw() +       
        theme(axis.text=element_text(size=20),axis.title=element_text(size=20),plot.title = element_text(color="black", size=20, face="bold")) +
        theme(legend.title = element_blank(),legend.text = element_text(size=20),legend.position = c(0.2,0.9),legend.direction = "horizontal")
    ggsave(paste0(cell[i],"_dens_distTADs.png"),width=10,height=10)
    ggsave(paste0(cell[i],"_dens_distTADs.pdf"),width=10,height=10)
    p <- c(p,wilcox.test(abs(dat$V8),abs(dat_nonsip$V8))$p.value)


    plotd <- data.table(plotd)
      ggplot(plotd, aes(x=label, y = abs(kb), fill=label)) +
        geom_violin(trim=TRUE) + scale_fill_brewer(palette="Paired") + ylab("Distance from SIP/non-SIP to TADs (KB)") + xlab("") + ylim(0,2000)+
        geom_point(data=plotd[, .("median"=median(abs(kb))), by="label"], aes(y=median), show.legend = F,size = 5)+ggtitle(cell[i]) +
        theme_bw() +       
        theme(axis.text=element_text(size=20),axis.title=element_text(size=20),plot.title = element_text(color="black", size=20, face="bold")) +
        theme(legend.title = element_blank(),legend.text = element_text(size=20))
    ggsave(paste0(cell[i],"_dens_distTADs_violon.png"),width=10,height=10)
    ggsave(paste0(cell[i],"_dens_distTADs_violon.pdf"),width=10,height=10)
}
p.res <- cbind(cell,format(p,scientific = T,digits=3))
#      cell               
# [1,] "Ery"    "2.19e-01"
# [2,] "MacMon" "1.18e-02"
# [3,] "Neu"    "2.85e-03"
# [4,] "MK"     "4.49e-02"
# [5,] "nCD4"   "5.04e-01"

###################
# 2. Gene density #
###################
setwd("/proj/yunligrp/users/jwen/SIP/revise/2_3")
library(ggplot2)
cell <- c("Ery","MacMon","Neu","MK","nCD4")
p <- NULL
for(i in 1:length(cell)){
    cat(cell[i],"\n")
    dat <- read.table(paste0(cell[i],"_SIP.geneDensity"),header = F)
    dat_nonsip <- read.table(paste0(cell[i],"_nonSIP.geneDensity"),header = F)
    dat$label="SIP"
    dat_nonsip$label <- "non-SIP"
    plotd <- rbind(dat,dat_nonsip)
    plotd <- data.table(plotd)
    ggplot(plotd, aes(x=label, y = V2, fill=label)) +
        geom_violin() + scale_fill_brewer(palette="Paired") + ylab("Gene density") + xlab("") +
        geom_point(data=plotd[, .("median"=median(V2)), by="label"], aes(y=median), show.legend = F,size = 5)+ggtitle(cell[i]) +
        theme_bw() +       
        theme(axis.text=element_text(size=20),axis.title=element_text(size=20),plot.title = element_text(color="black", size=20, face="bold")) +
        theme(legend.title = element_blank(),legend.text = element_text(size=20))
    ggsave(paste0(cell[i],"_geneDensity_violin.png"),width=10,height=8)
    ggsave(paste0(cell[i],"_geneDensity_violin.pdf"),width=10,height=8)
    p <- c(p,wilcox.test(dat$V2,dat_nonsip$V2)$p.value)
}
p.res <- cbind(cell,format(p,scientific = T,digits=3))
#      cell               
# [1,] "Ery"    "8.33e-02"
# [2,] "MacMon" "1.76e-10"
# [3,] "Neu"    "9.58e-06"
# [4,] "MK"     "1.93e-13"
# [5,] "nCD4"   "1.54e-14"