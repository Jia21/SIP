#!/bin/bash

# K562_sip=/proj/yunligrp/users/lagler/SIP/data/K562_SIPs.txt
# Note that  /proj/yunligrp/users/quansun/sip/K562_combined_r1r2r3.10k.2.peaks_SIPs.txt are the same with Taylor's calling results /proj/yunligrp/users/lagler/SIP/data/K562_SIPs.txt for K562 #
# Quan's results
dir_quan=/proj/yunligrp/users/quansun/sip
awk '$26==1' ${dir_quan}/K562_combined_r1r2r3.10k.2.peaks_SIPs.txt| awk '{print $2,$9,$10,$1}' OFS="\t" | sort | uniq > K562_Sips.bed 
awk '$26==0' ${dir_quan}/K562_combined_r1r2r3.10k.2.peaks_SIPs.txt| awk '{print $2,$9,$10,$1}' OFS="\t" | sort | uniq > K562_nonSips.bed

awk '$26==1' ${dir_quan}/GM_HiChIP_H3K27ac_b1b2_10kb_10k_2_SIPs.txt| awk '{print $2,$9,$10,$1}' OFS="\t" | sort | uniq > GM12878_Sips.bed
awk '$26==0' ${dir_quan}/GM_HiChIP_H3K27ac_b1b2_10kb_10k_2_SIPs.txt| awk '{print $2,$9,$10,$1}' OFS="\t" | sort | uniq > GM12878_nonSips.bed

gunzip ENCFF871HJT.bed.gz
cut -f1-3 ENCFF871HJT.bed > K562_hg38.bed
/proj/yunligrp/bin/liftOver K562_hg38.bed /proj/yunligrp/users/jwen/Reference/hg38ToHg19.over.chain.gz K562_CTCF.bed K562.unlifted
module add bedtools
bedtools intersect -a K562_Sips.bed -b K562_CTCF.bed -wa -wb > K562_SIP_CTCF
bedtools intersect -a K562_nonSips.bed -b K562_CTCF.bed -wa -wb > K562_nonSIP_CTCF

gunzip ENCFF710VEH.bed.gz
cut -f1-3 ENCFF710VEH.bed > GEM12878_CTCF_hg19.bed
bedtools intersect -a GM12878_Sips.bed -b GEM12878_CTCF_hg19.bed -wa -wb > GM12878_SIP_CTCF
bedtools intersect -a GM12878_nonSips.bed -b GEM12878_CTCF_hg19.bed -wa -wb > GM12878_nonSIP_CTCF

# 
cut -f4 GM12878_SIP_CTCF |sort | uniq -c |wc -l
# 594/823
cut -f4 GM12878_nonSIP_CTCF |sort | uniq -c |wc -l
# 4835/10274
# cut -f5-7 GM12878_SIP_CTCF |sort | uniq |wc -l
# 998/43631

cut -f4 K562_SIP_CTCF |sort | uniq -c |wc -l
# 540/811
cut -f4 K562_nonSIP_CTCF |sort | uniq -c |wc -l
# 5071/11352
# cut -f5-7 K562_SIP_CTCF |sort | uniq |wc -l
# 839/52064

# In R,
# K562 #
fisher.test(matrix(c(540,811-540,5071,11352-5071),ncol=2,byrow=F),alternative="greater")
# OR:2.46
# p-value < 2.2e-16

# GM12878 #
fisher.test(matrix(c(594,823-594,4835,10274-4835),ncol=2,byrow=F),alternative="greater")
# 2.92
# p-value < 2.2e-16