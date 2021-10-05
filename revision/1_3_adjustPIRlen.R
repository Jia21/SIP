setwd("/proj/yunligrp/users/jwen/SIP/revise/1_3")

# /proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/fraser.rda
library(data.table)
library("questionr")
library(vroom)
load("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[[31]] <- "MacMon"
library(tidyverse)
sczGWAS <- vroom("/proj/yunligrp/users/jwen/GWAS_sumstat/sullivan-pardinas2018/CLOZUK_PGC2noclo.METAL.assoc.dosage.fix",col_names = T)
sczGWAS_sigVars <- sczGWAS %>% filter(P < 5e-8) %>% dplyr::select(SNP, CHR, BP)
sczGWAS_sigVars$start <- sczGWAS_sigVars$BP-1
sczGWAS_sigVars$end <- sczGWAS_sigVars$BP
sczGWAS_sigVars$Phenotype <- "SCZ"
sczGWAS_sigVars$ancestry <- "EUR"
gwas_ct <- sczGWAS_sigVars

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

    # Adjust total PIR length #
    # fraser_ct$oelen <- fraser_ct$oeEnd - fraser_ct$oeStart + 1
    # Adjust SIPlength #
    fraser_ct$oelen <- fraser_ct$baitEnd - fraser_ct$baitStart + 1

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
write.table(mat,"sczGWASannot_adjustSIPlen_tables.txt",quote = F,sep="\t",col.names =T,row.names =F)

mat <- read.table("sczGWASannot_adjustSIPlen_tables.txt",header = T)
plotdata <- data.frame(mat)
plotdata$CT <- factor(plotdata$CT, levels=c("Neu", "nCD4", "MK", "MacMon", "Ery"))
plotdata$OR  <- as.numeric(as.matrix(plotdata$OR))
plotdata$pvalue  <- as.numeric(as.matrix(plotdata$pvalue))
plotdata$lower  <- as.numeric(as.matrix(plotdata$lower))
plotdata$upper  <- as.numeric(as.matrix(plotdata$upper))
allsipOR <- ggplot(plotdata[plotdata$type=="all_SIPs",], aes(x = OR, y = CT,label = formatC(pvalue,format = "e", digits = 2))) + 
  geom_point(size = 3.5, color = "purple3") + geom_text(vjust = -1.5,  hjust = 0.6, size = 4) +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = upper, xmin = lower),
                 size = .5, height = .2, color = "gray50") + 
  # facet_wrap(~type, scales="free", ncol=1,
  #            labeller=labeller(type=c(all.SIP="SIP", spec.SIP="Specific SIP"))) +
  theme_bw()+ theme(panel.grid.minor = element_blank()) +
  ylab("Cell Type Group") + xlab("Odds Ratio for Variant Overlap") 
allsipOR
ggsave("sczgwas_adjustSIPlen_ORs_allSIP.png", width=6, height=4, unit="in")


## For blood traits GWAS ##
vdat <- fread("/proj/yunligrp/users/lagler/SIP/data/ST3_ST4_BCX_Vuckovic_et_al.csv", drop=7)
vdat$ancestry <- "EUR"
# match phenotypes to chen
vdat[Phenotype=="RBC#", "Phenotype"] <- "RBC"
vdat[Phenotype=="WBC#", "Phenotype"] <- "WBC"
vdat[Phenotype=="NEUT#", "Phenotype"] <- "NEU"
vdat[Phenotype=="MONO#", "Phenotype"] <- "MONO"
vdat[Phenotype=="LYMPH#", "Phenotype"] <- "LYM"
vdat[Phenotype=="BASO#", "Phenotype"] <- "BASO"
vdat[Phenotype=="EO#", "Phenotype"] <- "EOS"
vdat[Phenotype=="PLT#", "Phenotype"] <- "PLT"

vdat2 <- vdat[Phenotype %in% c("HCT", "HGB", "MCH", "MCHC", "RBC", "RDW",
                           "MONO", "WBC", "PLT", "MPV", "LYM", "NEU") &
                Chr %in% 1:22]
colnames(vdat2)[6] <- "pval"
vdat2$pval <- 10^(-vdat2$pval) # match chen
vdat2$data <- "Vuckovic"

# Chen data has 5 ancestries Get union of SNPs
# /proj/yunligrp/users/jwen/BCX_meta_anal/GWASknownVars/Chen_*_2020.txt
files <- list.files(path="/proj/yunligrp/users/jwen/BCX_meta_anal/GWASknownVars",pattern="^Chen.*\\.txt$")
files <- paste0("/proj/yunligrp/users/jwen/BCX_meta_anal/GWASknownVars/",files)
chen <- lapply(files, fread)
anc <- tstrsplit(files, "_")[[2]] %>% unlist()

# columns and column names not all the same 
chen[[1]] <- chen[[1]][, c(1:4,11,15)]
colnames(chen[[1]])[1] <- "Pheno"
chen[[2]] <- chen[[2]][, c(1:4,11,15)]
chen[[3]] <- chen[[3]][, c(1:4,11,15)]
chen[[4]] <- chen[[4]][, c(1:4,11,14)]
colnames(chen[[4]])[4] <- "POS_(hg19)"
chen[[5]] <- chen[[5]][, c(1:4,9,15)]
# add ancestry
for(i in 1:5){
  chen[[i]]$ancestry <- anc[i]
  chen[[i]]$`P-value` <- as.numeric(chen[[i]]$`P-value`)
}
# merge all CTs together
chenall <- Reduce(function(...) merge(..., all = TRUE), chen)
colnames(chenall) <- c("Phenotype", "VariantID", "Chr", "Position", "pval",
                       "rsID", "ancestry")
chenall$data <- "Chen"

# combine sets 
gwasLong <- rbind(vdat2, chenall)
gwas <- gwasLong
gwas$Chr <- as.numeric(gwas$Chr)

RBC <- c("HCT", "HGB", "MCH", "MCHC", "RBC", "RDW")

CT <- c("Ery","MK","MacMon","Neu","nCD4")
# mat[,1] <- CT
mat <- NULL
windowsize=0
for (i in 1:length(CT)){
    cat(CT[i],"\n")
    if(CT[i]=="Ery"){ phenotype <- RBC }
    if(CT[i]=="MacMon"){phenotype <- "WBC" }
    if(CT[i]=="MK"){phenotype <- c("PLT", "MPV") }
    if(CT[i]=="nCD4"){phenotype <- c("LYM", "WBC") }
    if(CT[i]=="Neu"){phenotype <- c("NEU", "WBC") }

    gwas_ct <- gwas[Phenotype %in% phenotype,]

    setDT(gwas_ct)[, rsStart := Position - windowsize]
    setDT(gwas_ct)[, rsEnd := Position + windowsize]
    
    setkey(gwas_ct, Chr, rsStart, rsEnd)
    
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
    
    # Adjust total PIR length #
    # fraser_ct$oelen <- fraser_ct$oeEnd - fraser_ct$oeStart + 1
    # Adjust SIPlength #
    fraser_ct$oelen <- fraser_ct$baitEnd - fraser_ct$baitStart + 1

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
fwrite(mat, "bloodGWASannotSIP_adjustlen_tables.txt", sep="\t", col.names=T)

plotdata <- data.frame(mat)
plotdata$CT <- factor(plotdata$CT, levels=c("Neu", "nCD4", "MK", "MacMon", "Ery"))
plotdata$OR  <- as.numeric(as.matrix(plotdata$OR))
plotdata$pvalue  <- as.numeric(as.matrix(plotdata$pvalue))
plotdata$lower  <- as.numeric(as.matrix(plotdata$lower))
plotdata$upper  <- as.numeric(as.matrix(plotdata$upper))
allsipOR <- ggplot(plotdata[plotdata$type=="all_SIPs",], aes(x = OR, y = CT,label = formatC(pvalue,format = "e", digits = 2))) + 
  geom_point(size = 3.5, color = "purple3") + geom_text(vjust = -1.5,  hjust = 0.6, size = 4) +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = upper, xmin = lower),
                 size = .5, height = .2, color = "gray50") + 
  # facet_wrap(~type, scales="free", ncol=1,
  #            labeller=labeller(type=c(all.SIP="SIP", spec.SIP="Specific SIP"))) +
  theme_bw()+ theme(panel.grid.minor = element_blank()) +
  ylab("Cell Type Group") + xlab("Odds Ratio for Variant Overlap") 
allsipOR
ggsave("bloodgwasSIP_ORs_allSIP_label.png", width=6, height=4, unit="in")


allsipOR <- ggplot(plotdata[plotdata$type=="spec.SIPs",], aes(x = OR, y = CT,label = formatC(pvalue,format = "e", digits = 2))) + 
  geom_point(size = 3.5, color = "purple3") + geom_text(vjust = -1.5,  hjust = 0.6, size = 4) +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = upper, xmin = lower),
                 size = .5, height = .2, color = "gray50") + 
  # facet_wrap(~type, scales="free", ncol=1,
  #            labeller=labeller(type=c(all.SIP="SIP", spec.SIP="Specific SIP"))) +
  theme_bw()+ theme(panel.grid.minor = element_blank()) +
  ylab("Cell Type Group") + xlab("Odds Ratio for Variant Overlap") 
allsipOR
ggsave("bloodgwasSIP_ORs_specSIP_label.png", width=6, height=4, unit="in")
