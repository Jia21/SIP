setwd("/proj/yunligrp/users/jwen/SIP/revise/1_3")

# /proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/fraser.rda
library(data.table)
library("questionr")
library(vroom)
load("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[[31]] <- "MacMon"

atac_ery <- fread("Ery_peaks.narrowPeak.gz", select=c(1:3))
atac_mon <- fread("Mono_peaks.narrowPeak.gz", select=c(1:3))
atac_ncd4 <- fread("CD4_peaks.narrowPeak.gz", select=c(1:3))
atac_mk <- fread("Mega_peaks.narrowPeak.gz", select=c(1:3))


CT <- c("Ery","MK","MacMon","Neu","nCD4")
# mat[,1] <- CT
mat <- NULL
windowsize=0
for (i in 1:length(CT)){
    cat(CT[i],"\n")
    setDT(gwas_ct)[, rsStart := BP - windowsize]
    setDT(gwas_ct)[, rsEnd := BP + windowsize]
    
    setkey(gwas_ct, CHR, rsStart, rsEnd)
    
    # SIP cumulative interaction score
    # /proj/yunligrp/users/lagler/SIP/data/
    sipScore <- fread(paste0("/proj/yunligrp/users/lagler/SIP/data/SIP_", CT[i], "_specific.txt"), drop=3)
    colnames(sipScore)[[4]] <- "Specific"
    setorder(sipScore, -intscore)
    sipScore$rank <- seq(1:nrow(sipScore)) 
    fraser_ct <- fraser[eval(as.name(CT[i])) >= 5,]

     # combine fraser with sip and filter out non signficant interactions
    fraser_ct <- fraser_ct[, c(1:10)]
    bait_uniq_all <- data.frame(unique(fraser_ct$baitID))

    sip <- merge(sipScore, 
                fraser_ct,
                by="baitID")
    sip$oeChr <- as.numeric(sip$oeChr)

    setkey(sip, oeChr, oeStart, oeEnd)
    overlaps <- foverlaps(sip, gwas_ct, nomatch = NULL)
    
    baitID_GWAS <- unique(sort(overlaps$baitID))

    bait_uniq_all$label <- 0
    bait_uniq_all$label[bait_uniq_all$unique.fraser_ct.baitID. %in% baitID_GWAS] <- 1

    fraser_ct$oelen <- fraser_ct$oeEnd - fraser_ct$oeStart + 1

    t <- fraser_ct %>% group_by(baitID) %>% summarise(pir_len=sum(oelen))

    sip_ss <- sipScore$baitID[sipScore$Specific == 1]
    sip_as <- sipScore$baitID[sipScore$SIP == 1]

    bait_uniq_all$anno_ss <- 0
    bait_uniq_all$anno_ss[bait_uniq_all$unique.fraser_ct.baitID. %in% sip_ss] <- 1

    bait_uniq_all$anno_as <- 0
    bait_uniq_all$anno_as[bait_uniq_all$unique.fraser_ct.baitID. %in% sip_as] <- 1

    dat <- full_join(bait_uniq_all,t,by=c("unique.fraser_ct.baitID." ="baitID"))

    mod1 <- glm(label ~ pir_len + anno_as, data = dat, family = "binomial")
    or1 <- formatC(exp(mod1$coefficients["anno_as"]),format= "e",digits=2)
    p1 <- formatC(summary(mod1)$coefficients["anno_as","Pr(>|z|)"],format= "e",digits=2)
    lower1 <- formatC(odds.ratio(mod1, level=0.95)["anno_as",2],format= "e",digits=2)
    upper1 <- formatC(odds.ratio(mod1, level=0.95)["anno_as",3],format= "e",digits=2)

    mod2 <- glm(label ~ pir_len + anno_ss, data = dat, family = "binomial")
    or2 <- formatC(exp(mod2$coefficients["anno_ss"]),format= "e",digits=2)
    p2 <- formatC(summary(mod2)$coefficients["anno_ss","Pr(>|z|)"],format= "e",digits=2)
    lower2 <- formatC(odds.ratio(mod2, level=0.95)["anno_ss",2],format= "e",digits=2)
    upper2 <- formatC(odds.ratio(mod2, level=0.95)["anno_ss",3],format= "e",digits=2)

    mat_tmp <- rbind(c(CT[i],"all_SIPs",or1,p1,lower1,upper1),c(CT[i],"spec.SIPs",or2,p2,lower2,upper2))
    colnames(mat_tmp) <- c("CT","type","OR","pvalue","lower","upper")
    mat <- rbind(mat,mat_tmp)
}
write.table(mat,"sczGWASannot_adjustPIRlen_tables.txt",quote = F,sep="\t",col.names =T,row.names =F)
