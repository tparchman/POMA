#!/bin/bash

##########################################################################################
## POMA1 cleaning lane one
##########################################################################################

/working/parchman/tapioca/src/tap_contam_analysis --db /working/parchman/contaminants/illumina_oligos --pct 20 /working/parchman/POMA/POMA_S110_L002_R1_001.fastq > POMA.readstofilter.ill.txt 

echo "Illumina filtering done for lane 1"

/working/parchman/tapioca/src/tap_contam_analysis --db /working/parchman/contaminants/phix174 --pct 80 /working/parchman/POMA/POMA_S110_L002_R1_001.fastq > POMA.readstofilter.phix.txt 

echo "PhiX filtering done for lane 1"


/working/parchman/tapioca/src/tap_contam_analysis --db /working/parchman/contaminants/ecoli-k-12 --pct 80 /working/parchman/POMA/POMA_S110_L002_R1_001.fastq > POMA.readstofilter.ecoli.txt

echo "ecoli filtering done for lane 1"


cat /working/parchman/POMA/POMA_S110_L002_R1_001.fastq | fqu_cull -r POMA.readstofilter.ill.txt POMA.readstofilter.phix.txt POMA.readstofilter.ecoli.txt > POMA.clean.fastq

echo "Clean copy of lane 1 done"


