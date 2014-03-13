#!/bin/bash

mapping="bwa"
min_trimming=15
extr_method="c"
input="chr20_dwgsim_25cov_3000ins_incl35_sureselect_noerr_xsq_valid-1-1.bam"
out=$(echo ${input} | cut -d '.' -f1)
out=${out}_15minclip2sup
samtools_dir=""

java -Xmx4g -jar ./PotentialSplitMobileExtractor.jar -extr "$extr_method" -in "$input" -minh "$min_trimming" -out "$out" -tool "$mapping"

bwa aln -n 2 -l 1000 -t 4 -R 80 -k 0  -f "${out}_potentialsplitmob.sai" -c 52mobile_bwa_color "${out}_potentialsplitmob.fastq"
bwa samse -f "${out}_potentialsplitmob.sam" 52mobile_bwa_color "${out}_potentialsplitmob.sai" "${out}_potentialsplitmob.fastq"
 
java -Xmx4g -jar ./ExtractSplitAnchors.jar -mobbam "${out}_potentialsplitmob.sam" -biobam "$input" -tool "$mapping" -out "$out"

eval "${samtools_dir}samtools sort ${out}_splitanchors.bam ${out}_splitanchors_csort"
eval "${samtools_dir}samtools index ${out}_splitanchors_csort.bam"

java -Xmx8g -jar ./SplitAnchorClusterer.jar -in "${out}_splitanchors_csort.bam" -index "${out}_splitanchors_csort.bam.bai" -out "${out}" -minsup 2
