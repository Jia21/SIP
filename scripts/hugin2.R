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

setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")

# pcHiC data formatted for intra-chromosomal interactions and distance < 2MB
# /proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/
load("fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[[31]] <- "MacMon"

# BLUEPRINT expression data (mapps baitIDs to ENSEMBL ids)
# /proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/
expr <- fread("GeneExpressionMatrix_combine_Mac_Mon.txt") # Jia's file
#----------------------------------------------------------------------------------------------#
# EXAMPLE 1: MK PLT
#----------------------------------------------------------------------------------------------#
# sip subnetwork data
# "/proj/yunligrp/users/lagler/SIP/data/"
mk <- read.xlsx2("SIPsubnetworks.xlsx", sheetName = "MK_PLT")%>%
  subset(., baitID=="492452") %>% merge(., fraser[,c("baitID", "oeID", "MK", "dist")],
            by=c("baitID", "oeID"), all.x=T, all.y=F) %>% data.table()
mk$oeStart <- as.numeric(mk$oeStart)
mk$oeEnd <- as.numeric(mk$oeEnd)
mk$baitStart <- as.numeric(mk$baitStart)
mk$baitEnd <- as.numeric(mk$baitEnd)
mk$MK <- as.numeric(mk$MK)
mk$POS_.hg19. <- as.numeric(mk$POS_.hg19.)

# interaction level data
mkU <- unique(mk, by="oeID")
setDT(mkU)[, xstart := (baitStart + baitEnd)/2]
setDT(mkU)[, xend := (oeStart + oeEnd)/2]
setorder(mkU, dist)
mkU$curvature <- c(-.15, -.1, -.05)

# gene expression data
gene <- data.table("name"=c("CLCN2", "POLR2H", "CHRD", "THPO", "EPHB3"),
                   "start"=c(184079439, 184079506, 184097861, 184095932, 184279572),
                   "end"=c(184063973, 184086384, 184108524, 184089723, 184300197),
                   "ensemblID"=c("ENSG00000114859", "ENSG00000163882",
                                 "ENSG00000090539", "ENSG00000090534", "ENSG00000182580"),
                   "expr"=c(0.6218563,25.56251,0.0218635,0.001710853, 0.0005911997))

# atac-seq data
# /proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/Mega_peaks.narrowPeak.gz
atacmk <- fread("../../ATACseq/Mega_peaks.narrowPeak.gz", select = c(1:3))
colnames(atacmk) <- c("chr", "start", "end")
atacmk <- atacmk[chr == "chr3"]
setkey(mkU, oeStart, oeEnd)
atacmk_ovlp <- foverlaps(atacmk, mkU, by.x=c("start", "end"), nomatch=NULL) %>%
  unique(., by=c("start", "end"))


#----------------------------------------------------------------------------------------------#

# chicago score interactions
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
# #rsid <- ggplot(data=mk) + 
#   geom_segment(aes(x=POS_.hg19., xend=POS_.hg19., y=c(0:7), yend=c(1:8)))+
#   geom_label(aes(x=POS_.hg19., y=c(0.5:7.5), label=rsID)) + 
#   xlim(min(mk$oeStart)-7500, max(mk$baitEnd)+20000) +
#   theme_void() + theme(panel.background=element_rect(colour="black", fill=NA),
#                        axis.title.y = element_text(angle=90)) + 
#   ylab("Relevant Variants")
  rsid <- 
  ggplot(data=mk) + 
    geom_segment(aes(x=POS_.hg19., xend=POS_.hg19., y=0, yend=.25))+
    xlim(min(mk$oeStart)-7500, max(mk$baitEnd)+20000) +
    theme_void() + theme(panel.background=element_rect(colour="black", fill=NA),
                         axis.title.y = element_text(angle=90, size=7),
                         plot.margin = unit(c(11, 5.5, 0, 5.5), "points")) + 
    ylab("SNPs")
  
# gene expression
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
ggsave("../figures/HuginMKPLT.png", height=6, width=12)

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
# EXAMPLE 2: NCD4 LYM
#----------------------------------------------------------------------------------------------#
# sip subnetwork data
# setwd("/proj/yunligrp/users/lagler/SIP/data/")
nc <- read.xlsx2("SIPsubnetworks.xlsx", sheetName = "NCD4_LYM")%>%
  subset(., baitID=="139778") %>% merge(., fraser[,c("baitID", "oeID", "nCD4", "dist")],
                                        by=c("baitID", "oeID"), all.x=T, all.y=F) %>% data.table()
nc$oeStart <- as.numeric(nc$oeStart)
nc$oeEnd <- as.numeric(nc$oeEnd)
nc$baitStart <- as.numeric(nc$baitStart)
nc$baitEnd <- as.numeric(nc$baitEnd)
nc$nc <- as.numeric(nc$nc)
nc$POS_.hg19. <- as.numeric(nc$POS_.hg19.)

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
atacnc <- fread("../../ATACseq/CD4_peaks.narrowPeak.gz", select = c(1:3))
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
  #geom_label(aes(x=POS_.hg19., y=c(0.5, 0.5, 0.5), label=rsID)) + 
  xlim(min(nc$oeStart), max(nc$oeEnd)) + 
  #ylim(-.5,1.5)+
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
ggsave("../figures/HuginNCD4LYM.png", height=5, width=12)


#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
# both plots
ggarrange(mkplot, ncd4Plot, nrow=2,
        labels=letters[1:2], vjust=1, hjust=0, 
        heights = c(.525, .475), align="hv")
ggsave("../figures/Figure4_subnetworks2.png", height=8, width=8)
ggsave("../figures/Figure4_subnetworks2.pdf", height=8, width=8)

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
