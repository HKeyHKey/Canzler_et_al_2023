#!/bin/bash

sample=$1
adapter=NNNNTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG
cutadapt -q 20 -a $adapter --discard-untrimmed -m 18 -M 30 $sample'.fastq.gz' | cutadapt -u 4 - > $sample'_cutadapt.fastq'

