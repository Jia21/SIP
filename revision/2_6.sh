#!/bin/bash

k562_dir="/proj/yunligrp/users/jwen/allelic_imbalance/K562_genotypes/ENCFF752OAX.vcf.gz"
cd /proj/yunligrp/users/jwen/SIP/revise/2_6

module add vcftools
vcftools --gzvcf ${k562_dir}  --remove-indels --recode --recode-INFO-all --out ENCFF752OAX_SNPs_only
vcftools --vcf ENCFF752OAX_SNPs_only.recode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out ENCFF752OAX_biallele_SNPs_only
grep -v "#" ENCFF752OAX_biallele_SNPs_only.recode.vcf | cut -f1-2,10 | awk '{split($3,a,":"); {print $1,$2-1,$2,a[1]}}' OFS="\t" > K562_geno 

K562_sip=/proj/yunligrp/users/lagler/SIP/data/K562_SIPs.txt

awk '$26==1' $K562_sip |  cut -f11-13,22 |sort | uniq > K562_SIP_PIR
awk '$26==0' $K562_sip |  cut -f11-13,22 |sort | uniq > K562_nonSIP_PIR
module add bedtools
bedtools intersect -a K562_SIP_PIR -b K562_geno -wa -wb > K562_SIP_PIR_geno
bedtools intersect -a K562_nonSIP_PIR -b K562_geno -wa -wb > K562_nonSIP_PIR_geno

setwd("/proj/yunligrp/users/jwen/SIP/revise/2_6")
sip_pir_geno <- read.table("K562_SIP_PIR_geno",header = F)
non_sip_pr_geno <- read.table("K562_nonSIP_PIR_geno",header = F)
library(data.table)
library(tidyverse)

sip_pir_geno$V8[sip_pir_geno$V8=="0|1"] ="0/1"
sip_pir_geno$V8[sip_pir_geno$V8=="1|0"] ="0/1"
dat <- sip_pir_geno %>% select(V5,V6,V7,V8) %>% distinct() %>% group_by(V8) %>% 
        count(n=n()) %>% filter(V8!=1)
dat$label <- "SIP PIR"

non_sip_pr_geno$V8[non_sip_pr_geno$V8=="0|1"] ="0/1"
non_sip_pr_geno$V8[non_sip_pr_geno$V8=="1|0"] ="0/1"
dat_nosip_pir <- non_sip_pr_geno %>% select(V5,V6,V7,V8) %>% distinct() %>% group_by(V8) %>% 
        count(n=n()) %>% filter(V8!=1)
dat_nosip_pir$label <- "non-SIP PIR"
library(ggplot2)

plot <- rbind(dat,dat_nosip_pir)
# V8         n     nn label           
# 0/1   113309 113309 SIP PIR    
# 1|1   122847 122847 SIP PIR    
# 0/1   261018 261018 non-SIP PIR
# 1|1   278793 278793 non-SIP PIR
p <- ggplot(plot,aes(x=V8,y=n,fill=label)) +
     geom_bar(stat="identity", position=position_dodge()) +
     scale_fill_brewer(palette="Paired", direction=1) + 
     ylab("Count of genotypes") + xlab("Genotype") + 
     scale_y_continuous(breaks = c(0, 1e05, 2e05),labels = c("0",expression(1^5),expression(2^5))) +
     theme_bw()+theme(legend.title = element_blank()) +
     theme(axis.text=element_text(size=20),axis.title=element_text(size=20)) +
     theme(legend.title = element_blank(),legend.text = element_text(size=20))
ggsave("K562_geno.pdf",width=12,height=10)
ggsave("K562_geno.png",width=12,height=10)

## Extract 1KG EUR common variants ##
dir=/proj/yunligrp/projects/1000G/phase_3/freq/EUR

for i in {1..22}; do
        echo ${i}
        zcat ${dir}/chr${i}.frq.gz | sed 's/ \+/\t/g'  | sed 's/^\t\+//g' | awk '$5>0.01' |cut -f1-2  |\
        sed 1d | awk '{split($2,a,":"); {print "chr"$1,a[2]-1,a[2],$1":"a[2]":"a[3]":"a[4]}}' OFS="\t" > \
                /proj/yunligrp/users/jwen/SIP/revise/2_6/1KG_EUR/chr${i}.commonvars.frq
done

cat chr*.commonvars.frq > chrs.commonvars.frq
awk '$4=="normal"' /proj/yunligrp/users/jwen/allelic_imbalance/copy_variation/ENCFF000LPC.bed  |  cut -f1-3 | awk '$1!="chrX"' > K562_normal_regions

module add bedtools
bedtools intersect -wao -a /proj/yunligrp/users/jwen/SIP/revise/2_6/1KG_EUR/chrs.commonvars.frq -b K562_normal_regions | awk '$5!="."' > allvars_KG_EUR_normal_region
# 6449616

bedtools intersect -wao -a K562_geno -b K562_normal_regions | awk '$5!="."' > K562_geno_normal_region

bedtools intersect -a K562_SIP_PIR -b K562_geno_normal_region -wa -wb > K562_SIP_PIR_geno_normal_region
cut -f 5-8 K562_SIP_PIR_geno_normal_region | sort | uniq -c|wc -l
# 172,319 geno in sip pir region and in k562 normal region
cut -f5-8 K562_SIP_PIR_geno_normal_region | sort | uniq |cut -f4 |sort | uniq -c 
#   42283 0|1
#    1755 0/1
#   42300 1|0
#   85981 1|1
# 580,120 - 172,319 = 407,801 0|0
bedtools intersect -a K562_nonSIP_PIR -b K562_geno_normal_region -wa -wb > K562_nonSIP_PIR_geno_normal_region
cut -f 5-8 K562_nonSIP_PIR_geno_normal_region | sort | uniq -c|wc -l
# 397,513 geno in nonsip pir region and in k562 normal region
cut -f5-8 K562_nonSIP_PIR_geno_normal_region | sort | uniq |cut -f4 |sort | uniq -c 
#   99197 0|1
#    3874 0/1
#  102608 1|0
#  191834 1|1
# 1,327,890 - 397,513 = 930,377 0|0

bedtools intersect -wao -a allvars_KG_EUR_normal_region -b K562_SIP_PIR |awk '$9!="."' > K562_SIP_PIR_allvars_KG_EUR_normal_region
cut -f 1-4 K562_SIP_PIR_allvars_KG_EUR_normal_region | sort | uniq -c|wc -l
# 1kg 580,120 common vars in sip pir region, k562 nomral region 
bedtools intersect -wao -a allvars_KG_EUR_normal_region -b K562_nonSIP_PIR |awk '$9!="."' > K562_nonSIP_PIR_allvars_KG_EUR_normal_region
cut -f 1-4 K562_nonSIP_PIR_allvars_KG_EUR_normal_region | sort | uniq -c|wc -l
# 1kg 1,327,890 common vars in sip pir region, k562 nomral region 

library(ggplot2)

# V8         n     nn label  
# 0|0   407801  407801  SIP PIR
# 0/1   86338 86338 SIP PIR    
# 1|1   85981 85981 SIP PIR
# 0|0   930377  930377  non-SIP PIR
# 0/1   205679 205679 non-SIP PIR
# 1|1   191834 191834 non-SIP PIR
plot <- cbind(c("0/0","0/1","1/1","0/0","0/1","1/1"),c(407801/580120,86338/580120,85981/580120,930377/1327890,205679/1327890,191834/1327890),c(rep("SIP PIR",3),rep("non-SIP PIR",3)))
plot <- data.frame(plot)
p <- ggplot(plot,aes(x=X1,y=X2,fill=X3)) +
     geom_bar(stat="identity", position=position_dodge()) +
     scale_fill_brewer(palette="Paired", direction=1) + 
     ylab("Ratio of genotypes") + xlab("") + 
     theme_bw()+theme(legend.title = element_blank()) +
     theme(axis.text=element_text(size=20),axis.title=element_text(size=20)) +
     theme(legend.title = element_blank(),legend.text = element_text(size=20))
ggsave("./2_6/K562_geno_ratio.pdf",width=12,height=10)
ggsave("./2_6/K562_geno_ratio.png",width=12,height=10)