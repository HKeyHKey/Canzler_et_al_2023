# Canzler_et_al_2023
Data and scripts for the analysis of small RNA-Seq data presented in [Canzler et al., 2023](https://pubmed.ncbi.nlm.nih.gov/39441382/).

## Preparation of hisat2 index with extended rat pre-miRNA hairpin sequences: ##

``./Script_build_index.sh Extended_Improved_rat_pre-miRNA``

## Adapter trimming: ##

``for sample in `cat sample_list`;do ./Script_trimming_adapter.sh $sample;done``

## Mapping trimmed reads on extended pre-miRNA sequences: ##

``for sample in `cat sample_list`;do ./Script_hisat2.sh $sample Extended_Improved_rat_pre-miRNA;done``

## Counting miRNA isoform abundance: ##

``for sample in `cat sample_list`;do ./Script_count_miRNA.sh $sample;done``
