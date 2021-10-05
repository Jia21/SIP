setwd("/proj/yunligrp/users/jwen/SIP/revise/2_8")
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(kableExtra)
#----------------------------------------------------------------------------------------------#
# find overlaps with ATAC-seq by CT
#----------------------------------------------------------------------------------------------#
load("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[31] <- "MacMon"
allSIPgeneInfo <- read.table("/proj/yunligrp/users/lagler/SIP/data/allSIPs_geneInfo.txt",fill=T,sep="\t",header =T)

atac_ery <- fread("/proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/Ery_peaks.narrowPeak.gz", select=c(1:3))
atac_mon <- fread("/proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/Mono_peaks.narrowPeak.gz", select=c(1:3))
atac_ncd4 <- fread("/proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/CD4_peaks.narrowPeak.gz", select=c(1:3))
atac_mk <- fread("/proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/Mega_peaks.narrowPeak.gz", select=c(1:3))

CT<-c("Ery","nCD4","MK","MacMon")
pl <- NULL
p <- NULL
for(i in 1:length(CT)){
  cat(CT[i],"\t")
  hic <- fraser %>% filter(!!as.name(CT[i]) >= 5) %>% data.table
  hic <- hic[,1:11] 
  sip <- fread(paste0("/proj/yunligrp/users/lagler/SIP/data/SIP_", CT[i], ".txt"))
  sipdata <- merge(hic, sip[sip$SIP==1,], by="baitID", allow.cartesian = T)

  setDT(sipdata)[, chr2 := paste0("chr", baitChr)]

  if(CT[i]=="Ery"){colnames(atac_ery) <- c("chr", "start", "end"); seqdat <- atac_ery;}
  if(CT[i]=="MacMon"){colnames(atac_mon) <- c("chr", "start", "end"); seqdat <- atac_mon;}
  if(CT[i]=="nCD4"){colnames(atac_ncd4) <- c("chr", "start", "end"); seqdat <- atac_ncd4;}
  if(CT[i]=="MK"){colnames(atac_mk) <- c("chr", "start", "end"); seqdat <- atac_mk;}

  seqdat$size <- seqdat$end - seqdat$start + 1
  seqdat$row <- seq(1:nrow(seqdat))
  setkey(seqdat, chr, start, end)
  overlap <- foverlaps(sipdata, seqdat, by.x=c("chr2", "baitStart", "baitEnd"),
                       type="any", mult="first", which=T)

  sipdata$atacRow <- overlap # row number of ATAC-seq overlap for size merging
  sipdata$atac <- as.numeric(is.na(overlap)==F) # recode to 0/1 (false/true) overlap
  
  sip_atac <- unique(sipdata$baitID[sipdata$atac == 1])
  sip_noatac <- unique(sipdata$baitID[sipdata$atac == 0])

  sip_atac_gene <- allSIPgeneInfo %>% filter(baitID %in% sip_atac) %>% 
                    select(baitID,CT[i]) %>% mutate(label="SIPs olap ATAC peaks") %>% 
                    mutate(CT=CT[i])
  sip_noatac_gene <- allSIPgeneInfo %>% filter(baitID %in% sip_noatac) %>% 
                    select(baitID,CT[i]) %>% mutate(label="SIPs no olap ATAC peaks") %>% 
                    mutate(CT=CT[i]) 
  p <- c(p,wilcox.test(na.omit(sip_atac_gene)[,2],na.omit(sip_noatac_gene)[,2],alternative="greater")$p.value)
   dat <- rbind(sip_atac_gene,sip_noatac_gene) %>% na.omit(dat)
   colnames(dat)[2] <- "exp"
   pl <- rbind(pl,dat)
}

library(ggplot2)
pl <- data.table(pl)
    plotd <- ggplot(pl, aes(x=label, y = log(exp), fill=label)) +
          geom_violin() + scale_fill_brewer(palette="Paired") + ylab("log(Gene expression)") + 
          xlab("") +
          geom_point(data=pl[, .("median"=median(log(exp))), by=label], aes(y=median), show.legend = F,size = 3)+ 
          facet_wrap(.~CT) + theme_bw() + 
          theme(axis.text=element_text(size=20),axis.title=element_text(size=20)) +
          theme(legend.title = element_blank(),legend.text = element_text(size=20))+
          theme(strip.text.x = element_text(size = 15))
 ggsave(paste0("./","SIPs_exp_violin.png"),width=20,height=10)
 ggsave(paste0("./","SIPs_exp_violin.pdf"),width=20,height=10)