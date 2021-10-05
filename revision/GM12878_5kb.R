#----------------------------------------------------------------------------------------------#
# Call SIPs in K562 HiChIP data
# requires K562 MAPS HiChIP data and gencode promoters
# defines SIPS, outputs K562_SIPs.txt
# computes basic SIP properties and tests if SIPs are more likely to be enhancers than non-SIPs
#----------------------------------------------------------------------------------------------#
library(data.table)
library(ggplot2)
library(dplyr)

#----------------------------------------------------------------------------------------------#

# MAPS calls on K562 HiChIP data from Ming
# mapsDir <- "/proj/yunligrp/MingHuData/H3K27ac_HiChIP_K562_MAPS_interactions/"
mapsDir <- "/proj/yunligrp/users/minghu/081521_HiChIP_data_for_SIP_revision/"
k <- fread(paste0(mapsDir, "GM_HiChIP_H3K27ac_b1b2_5kb.5k.2.sig3Dinteractions.bedpe"))

# hg19 promoters from gencode
# gencodeDir <- "/proj/yunligrp/users/ycyang/Gencode/"
gencodeDir <- "/proj/yunligrp/users/quansun/annotation/"
gencode <- fread(paste0(gencodeDir, "gene_info.gencode.v28.hg19.txt"),
                 drop=c("Strand", "Gene_Length"))

# overlap HiCHiP anchor with promoter
setkey(k, "chr1", "start1", "end1")
setkey(gencode, "CHR", "Promoter_Start", "Promoter_End")
b1P <- foverlaps(k, gencode, nomatch = NULL) # keep if bait 1 is promoter
setkey(k, "chr2", "start2", "end2")
b2P <- foverlaps(k, gencode, nomatch = NULL) # keep if bait 2 is promoter

# reformat data to be in bait/PIR style
b1P <- rename(b1P, c(baitChr=chr1, baitStart=start1, baitEnd=end1,
                     oeChr=chr2, oeStart=start2, oeEnd=end2))
b2P <- rename(b2P, c(baitChr=chr2, baitStart=start2, baitEnd=end2,
                     oeChr=chr1, oeStart=start1, oeEnd=end1))
# combine both
kP <- rbind(b1P, b2P)
setDT(kP)[, baitID := paste0(baitChr, ":", baitStart, ":", baitEnd)]
setDT(kP)[, oeID := paste0(oeChr, ":", oeStart, ":", oeEnd)]

#----------------------------------------------------------------------------------------------#
# use -log10 FDR as SIP score
setDT(kP)[, log10fdr := -log10(fdr)]
setDT(kP)[, log10fdr := ifelse(log10fdr==Inf, 320, log10fdr)]

summary(kP[, .N, by="baitID"]$N)

# sum scores by bait to get cumulative interaction score
score <- kP[, .(intscore = sum(log10fdr)), by="baitID"]
# rank cumulative interaction scores
setorder(score, intscore)
score$rank <- seq(1:nrow(score))

# find inflection point (i.e. define SIPs)
# normalize rank and score
setDT(score)[, norm_rank := rank/max(rank)]
setDT(score)[, norm_score := intscore/max(intscore)]
setDT(score)[, rot1 := 1/sqrt(2)*norm_rank + 1/sqrt(2)*norm_score]
setDT(score)[, rot2 := -1/sqrt(2)*norm_rank + 1/sqrt(2)*norm_score]

RefPoint <- score[rot2==min(rot2), rank]
RefValue <- score[rot2==min(rot2), intscore] # any score > is a SIP
# add SIP indicator
setDT(score)[, SIP := ifelse(intscore > RefValue, 1, 0)]
NumSIP <- sum(score$SIP)

# add SIP indicator to data
sips <- merge(kP, score[,c("baitID", "intscore", "rank", "SIP")], by="baitID")
fwrite(sips, "GM_HiChIP_H3K27ac_b1b2_5kb_5k_2_SIPs.txt", sep="\t", col.names=T)

table(unique(sips, by="baitID")$SIP) # 811 SIPs, 11,352 non-SIPs
table(unique(sips, by=c("baitID", "Gene_ID"))$SIP) # 1,348 SIP genes (1,284 unique); 13,215 non-SIP genes
#----------------------------------------------------------------------------------------------#
# are SIPs more likely to be an enhancer for another promoter bait than non-SIPs?
baitU_sip <- unique(sips[SIP==1]$baitID) #unique SIP baits
baitU_nonsip <- unique(sips[SIP==0]$baitID) #unique non-SIP baits
pirU <- unique(sips$"oeID") #unique PIRs

n.sipE <- length(intersect(baitU_sip, pirU)) # number of SIPs that are also enhancers/PIRs
n.nonsipE <- length(intersect(baitU_nonsip, pirU)) #number of non-SIPs that are also enhancers/PIRs

x <- matrix(c(n.sipE, length(baitU_sip)-n.sipE,
              n.nonsipE, length(baitU_nonsip)-n.nonsipE), ncol = 2, byrow=T)
colnames(x) <- c("Enh", "no Enh")
rownames(x) <- c("SIP", "non-SIP")
x
chisq.test(x)
fisher.test(x)

#----------------------------------------------------------------------------------------------#
# hockey stick plot
ggplot(score, aes(x=rank, y=intscore)) + 
  geom_point(size=1) + 
  geom_vline(xintercept=RefPoint, color="blue", linetype="dashed") +
  xlab("Rank of Promoter Bait") + ylab("Cumulative Interaction Score") + 
  #scale_x_continuous(breaks=c(0, 5000, 10000))+
  ggtitle("GM12878 H3K27ac HiChIP SIPs (5kb)")+
  theme_bw()
ggsave("./figures/GM_HiChIP_H3K27ac_b1b2_5kb_5k_2_SIPs.png", width=6, height=5)

# median score per bait
medScore <- sips[, .(M = median(log10fdr)), by="baitID"] %>% 
  merge(., sips[,c("baitID", "SIP")], by="baitID", all=F)  %>%
  unique(., by="baitID")

summary(medScore[SIP==1, M])
summary(medScore[SIP==0, M])
wilcox.test(medScore[SIP==1, M], medScore[SIP==0, M])

ggplot(data=medScore, aes(x=SIP, y=log10(M), fill=as.factor(SIP))) + 
  geom_violin()+
  geom_point(data=medScore[, .("median"=median(M)), by=c("SIP")],
             aes(y=log10(median)), position =position_dodge(.9), show.legend = F)+
  scale_fill_brewer(palette="Paired", direction=1, labels = c("non-SIP", "SIP"))+ 
  ylab("Median MAPS Score per Bait (log10)") + 
  xlab("Bait Type") + 
  theme_bw()+theme(legend.title = element_blank())
ggsave("./figures/GM_HiChIP_H3K27ac_b1b2_5kb_5k_2_score.png", width=5, height=6)

# # of PIRs per bait
nPIR <- sips[, .N, by="baitID"] %>% 
  merge(., sips[,c("baitID", "SIP")], by="baitID", all=F) %>%
  unique(., by="baitID")

summary(nPIR[SIP==1, N])
summary(nPIR[SIP==0, N])
wilcox.test(nPIR[SIP==1, N], nPIR[SIP==0, N])

ggplot(data=nPIR, aes(x=SIP, y=log10(N), fill=as.factor(SIP))) + 
  geom_violin()+
  geom_point(data=nPIR[, .("median"=median(N)), by=c("SIP")],
             aes(y=log10(median)), position =position_dodge(.9), show.legend = F)+
  scale_fill_brewer(palette="Paired", direction=1, labels = c("non-SIP", "SIP"))+ 
  ylab("Number of Interactions per Bait (log10)") + 
  xlab("Bait Type") + 
  theme_bw()+theme(legend.title = element_blank())
ggsave("./figures/GM_HiChIP_H3K27ac_b1b2_5kb_5k_2_PIR.png", width=5, height=6)
#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#

