#!/usr/bin/bash

# Written by SJK 04. 29. 2020
# Purpose: This file is for GEO fastq files download

# 1. Prerequirement
# To run the file, require two files SRR_list.txt & SRR_getdata.sh
# Manually filled out wanted SRR IDs in SRR_list.txt
# Download sratoolkit.2.9.0-ubuntu64  

# 2. How to Run
# ./SRR_download.sh $


# 3. Actual Run

cat SRR_list.txt |xargs -n 1 ./SRR_getdata.sh &

