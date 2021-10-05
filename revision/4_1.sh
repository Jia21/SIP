#!/bin/bash

awk '$26==1' GM12878_rep_1_SIPs.txt |cut -f1  | sed 1d | sort | uniq > GM12878_rep_1.SIPs
awk '$26==1' GM12878_rep_2_SIPs.txt |cut -f1  | sed 1d | sort | uniq > GM12878_rep_2.SIPs
comm -12 <(sort GM12878_rep_1.SIPs) <(sort GM12878_rep_2.SIPs) | wc -l
# 491
wc -l GM12878_rep_*.SIPs
#  596 GM12878_rep_1.SIPs
#  515 GM12878_rep_2.SIPs

for rep in {1..3};do 
    echo $rep
    awk '$26==1' K562_rep_${rep}_SIPs.txt |cut -f1  | sed 1d | sort | uniq > K562_rep_${rep}.SIPs
done