#----------------------------------------------------------------------------------------------#
# compares SIPs in K562 to validated CRISPR data from 3 papers
# finds overlaps between promoters or PIRs and gRNAs/disrupted regions
# 2x2 tables of results by significance/SIP status

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr) # separate
library(readxl) # read_xlsx
library(Biobase) # get gene names from ensembl ids
library(biomaRt)

setwd("/Users/tmlagler/OneDrive/Lab/SIP/figurecode/data/")
#setwd("/proj/yunligrp/users/lagler/SIP/data/")

# SIP calls on MAPS K562 HiChIP data
sips <- fread("K562_SIPS.txt")

#----------------------------------------------------------------------------------------------#
# Fulco et al. 2019 ABC CRISPR
# Supplemental Table 3a
fulco <- read_xlsx("Fulco2019supplement.xlsx",
                   sheet="Supplementary Table 3a",
                   skip=1) %>% data.table()

# only select promoter elements and valid connections
# chrs 3, 12, and 19 (SIP data doesn't have chrX, remove here so counts are accurate)
fulco <- separate(fulco, "Element name", c("type", "position"), "\\|") 
fulco <- fulco[type=="promoter" & `Valid E-G connection`==T
               & chr != "chrX"]

table(fulco$Significant) # 39 significant out of 920

# overlap disrupted enhancer with promoter bait
setkey(sips, "baitChr", "baitStart", "baitEnd")
setkey(fulco, "chr", "start", "end")
sips_fulco <- foverlaps(sips, fulco, nomatch=NULL)
# match target gene and SIP gene, take unique by promoter bait
# in this case, one SIP corresponds two two genes
sips_fulco <- sips_fulco[Gene == Gene_NAME] %>% unique(., by=c("baitID", "Gene_NAME"))
sips_fulcoU <- sips_fulco[Gene == Gene_NAME] %>% unique(., by="baitID")

table(sips_fulcoU$SIP, sips_fulcoU$Significant)
sips_fulco[SIP==1 & Significant==T, Gene_NAME]

#----------------------------------------------------------------------------------------------#
# Morris et al. bioRxiv 2021
# Table S3E  STING-seq cis-regulatory element results, tested within the SCEPTRE framework. 
# The top gRNAs for 37 variants identified 33 target genes in cis.
morris <- read_xlsx("Morris_TableS3.xlsx",
                    sheet="Table S3E",
                    skip=1) %>% data.table()

# define significant vs not-significant
setDT(morris)[, significant := ifelse(is.na(`Significance tier`)==T, FALSE, TRUE)]
table(morris$significant)

# find position of SNPs
snpInfo <- read_xlsx("morrisSNPinfo.xlsx") %>% data.table()
morris <- merge(morris, snpInfo, by="SNP") %>% unique()

# overlap SNPs with PIRs (enhancers)
setkey(sips, "oeChr", "oeStart", "oeEnd")
setkey(morris, "chr", "startSNP", "endSNP")
sips_morris_pir <- foverlaps(sips, morris, nomatch=NULL)  
# match target gene and SIP gene, take unique by promoter/PIR 
sips_morris_pir <- sips_morris_pir[Gene==Gene_NAME] %>% unique(., by=c("oeID", "baitID"))

table(sips_morris_pir$SIP, sips_morris_pir$significant)
sips_morris_pir[SIP==1 & significant==T, "Gene_NAME"]

# overlap SNPs with SIPs
setkey(sips, "baitChr", "baitStart", "baitEnd")
sips_morris_bait <- foverlaps(sips, morris, nomatch=NULL)  
# match target gene and SIP gene, take unique by promoter bait
sips_morris_bait <- sips_morris_bait[Gene==Gene_NAME] %>% unique(., by=c("baitID", "Gene_NAME"))
sips_morris_baitU <- unique(sips_morris_bait, by="baitID")

table(sips_morris_baitU$SIP, sips_morris_baitU$significant)
sips_morris_bait[SIP==1 & significant==T, "Gene_NAME"]

#----------------------------------------------------------------------------------------------#
# Klann et al. bioRxiv 2021
# Table S6 K562 validation screen gRNA info and raw counts and DESeq2 results
klann <- fread("supplementary_table_6_validation_screen_k562_sgrna_deseq2_results_hg19.csv.gz")
klann <- klann[chr!=NaN & chr!="chrX"]
# define significance
setDT(klann)[, significant := ifelse(padj < 0.05, TRUE, FALSE)]
table(klann$significant)

# overlap gRNA with SIP PIR
setkey(sips, "oeChr", "oeStart", "oeEnd")
setkey(klann, "chr", "start", "end")
sips_klann_pir <- foverlaps(sips, klann, nomatch=NULL)
# unique by promoter/PIR
sips_klann_pir <- unique(sips_klann_pir, by=c("baitID", "oeID"))

table(sips_klann_pir$SIP, sips_klann_pir$significant)
uniqueN(sips_klann_pir[SIP==1 & significant==T], by="baitID")
chisq.test(table(sips_klann_pir$SIP, sips_klann_pir$significant))

# overlap gRNA with promoter bait
setkey(sips, "baitChr", "baitStart", "baitEnd")
sips_klann_bait <- foverlaps(sips, klann, nomatch=NULL) 
# unique by promoter bait (i.e., at least one overlap)
sips_klann_bait <- unique(sips_klann_bait, by=c("baitID", "Gene_NAME"))
sips_klann_baitU <- unique(sips_klann_bait, by="baitID")

table(sips_klann_baitU$SIP, sips_klann_baitU$significant)
sips_klann_bait[SIP==1 & significant==T, Gene_NAME]


#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#

