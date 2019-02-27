#!/bin/bash

#####################################################
# 15 Feb 2019 by SJK
# This is file for report QC for single cell RNA seq
# Below tools was used to create QC metrics file
# Drop-seq pipeline version 2.0.0
# samtools version 1.9 
####################################################

# How to use
#./scRNA_seq_QC_metrics.sh ctx_12 > QC_Report.txt

ARG1=$1
##################
# 1. Flagstats
##################
echo "1st Flag stats of BAM file"
samtools flagstat my_clean_$ARG1'.bam'

#########################
# 2. MAPQ quality stats
#########################
# This stats can identify how many uniquely mapped reads
echo "2nd MAPQ stats"
samtools view my_clean_$ARG1'.bam' |awk '{h[$5]++};END { for(k in h) print k, h[k]}' > $ARG1'MAPQ.txt'
cat $ARG1'MAPQ.txt'

##################################################
# 3. Extract the proportion of each region mapping
##################################################
# for percentage calculation
calc() { awk "BEGIN{ printf \"%.3f\n\", $* }"; }
# The number of uniquely mapped reads
Uniq_R=`awk 'NR == 1 {print $2}' $ARG1'MAPQ.txt'`

# In STAR mapping, MAPQ 255 is uniquely mapped read  
# Calculate the percentage of where the reads mapped
echo "3rd. mapping regions of uniquely mapped reads"
UTR=`samtools view my_clean_$ARG1'.bam' |awk '$5 == 255' | grep -c "XF:Z:UTR"`
UTRP=`calc $UTR / $Uniq_R`
echo "UTR: " $UTR "(" $UTRP "%)" 
CODING=`samtools view my_clean_$ARG1'.bam' |awk '$5 == 255' | grep -c "XF:Z:CODING"`
CODINGP=`calc $CODING / $Uniq_R`
echo "CODING: " $CODING "(" $CODINGP "%)"
INTRONIC=`samtools view my_clean_$ARG1'.bam' |awk '$5 == 255' | grep -c "XF:Z:INTRONIC"`
INTRONICP=`calc $INTRONIC / $Uniq_R`
echo "INTRONIC: " $INTRONIC "(" $INTRONICP "%)"
INTERGENIC=`samtools view my_clean_$ARG1'.bam' |awk '$5 == 255' | grep -c "XF:Z:INTERGENIC"`
INTERGENICP=`calc $INTERGENIC / $Uniq_R`
echo "INTERGENIC: " $INTERGENIC "(" $INTERGENICP "%)"

#####################
# 4. PCR calculation
#####################
# This calculation is done before filtering nGene or nUMI
# median of PCR is returned for each sample
echo "4th. PCR calculation"
for file in *.summary.txt
do
	PCR=`awk 'NR > 7 {print}' $file |awk '{print($2/$3)}'|sort -n| awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'`
	echo $file
	echo "PCR exon only: " $PCR
done

# Remove unnecssary MAPQ file
rm $ARG1'MAPQ.txt'








