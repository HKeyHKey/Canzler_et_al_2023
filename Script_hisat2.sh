#!/bin/bash

sample=$1
index=$2
hisat2 -k 1 --no-unal --norc -x $index -U $sample'_cutadapt.fastq' -S Mapping_$sample'_trimmed_0.sam' --un Unmapping_$sample'_trimmed_0.fastq'
for step in `seq 1 5` # iterative trimming on the 3´ end (to detect 3´-tailed isoforms)
do previous=`echo $step"-1" | bc`
   hisat2 -3 $step -k 1 --no-unal --norc -x $index -U Unmapping_$sample'_trimmed_'$previous'.fastq' -S Mapping_$sample'_trimmed_'$step'.sam' --un Unmapping_$sample'_trimmed_'$step'.fastq'
done
