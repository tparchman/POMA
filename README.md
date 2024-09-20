
## Organization and Workflow for POMA (8/24) 
Organizational notes and code for one sequencing sets:
- 1 plus plates in the first freezer door have DNAs in order as loaded into 96 well plates. See also `plate_map.csv` 

## Sample organization
- Full information on DNAs for each individual sampled across natural distribution can be found in `GitHub/RFseq_projects/POMA/plate_map.csv`. This file also has the updated plate maps with specified IDs.

- barcode key file corresponds with POMA:`POMA_barcode_key.csv`


## Notes on library preparation


### R/L and PCR for plates 1-4. 

- Master mix in `PM_skipper_RFseq_mastermixcockatils.xlsx`.
- Restriction done on 5/2/24
- Ligation complete on 5/3/24
- PCR done on: 5/5/24



10 ul of each PCR product into final library. Tubes in door of freezer labelled **POMA**.


## Data analysis: contaminant cleaning, barcode parsing, data storage, directory organization, and initial analyses.

We generated one lane of S1 chemistry NovaSeq data at UTGSAF in July of 2024. 


## This file contains code and notes for
1) cleaning contaminants using tapioca
2) parsing barcodes
3) splitting fastqs 
4) de novo assembly
5) reference based assembly
6) calling variants
7) filtering
8) entropy for genotype probabilities.

## 1. Cleaning contaminants

Being executed on ponderosa using tapioca pipeline. Commands in two bash scripts (cleaning_bash_CADE_FRVE.sh), executed as below (4/12/24). This work being carried out at:

    /working/parchman/POMA/

Decompress fastq files:

    $ nohup gunzip POMA_S110_L002_R1_001.fastq.gz &
Number of reads **before** cleaning:

    $ nohup grep -c "^@" POMA_S110_L002_R1_001.fastq > POMA.txt &
    ## raw reads: 

To run cleaning_bash* tapioca wrapper, exit conda environment, load modules, and run bash scripts.

    $ module load fqutils/0.4.1
    $ module load bowtie2/2.2.5
    
    $ bash cleaning_bash_POMA.sh &


Raw data (POMA_S110_L002_R1_001.fastq) will stay stored in: /backups/rawdata_to_backup/POMA/

Number of reads **after** cleaning:

    $ nohup grep -c "^@" OVCA.clean.fastq > OVCA_clean_reads.txt &
    # number of clean reads : 


####################################################################################
## 2. Barcode parsing:
####################################################################################

Be sure to deactivate conda environment before running the below steps. Barcode keyfile is `/working/parchman/POMA/POMA_barcode_key.csv`
`
Parsing POMA library:

    $ nohup perl parse_barcodes768.pl POMA_barcode_key.csv POMA.clean.fastq A00 &>/dev/null &

looking at report:

    $ less parsereport_POMA.clean.fastq
    Good mids count: 594526790
    Bad mids count: 50658963
    Number of seqs with potential MSE adapter in seq: 227974
    Seqs that were too short after removing MSE and beyond: 125

####################################################################################
## 3. splitting fastqs
####################################################################################

Doing this in `/working/parchman/POMA/splitfastqs_POMA/` 

Make ids file

    $ cut -f 3 -d "," POMA_barcode_key.csv | grep "[A-Z]" > POMA_ids_noheader.txt

Split fastqs by individual

    $ nohup perl splitFastq_universal_regex.pl POMA_ids_noheader.txt parsed_POMA.clean.fastq &>/dev/null &



# DONE TO HERE &&&&&&&&&&


Zip the parsed*fastq files for now, but delete once patterns and qc are verified.

