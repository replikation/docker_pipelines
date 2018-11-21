#!/usr/bin/env bash

## Docker ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found. Please run the installation skript, Aborting.${NC}"; exit 1; }
  echo "Docker identified"

## OPTIONS ##
  WORKDIRPATH=$(pwd) # for docker mountpoint (-v)
  WORKDIRNAME=${PWD##*/} # for docker mountpoint (-v)
  # SCRIPTLOCATION=$(readlink -f "$0") # Script location
  # SCRIPTPATH=$(dirname "$SCRIPTLOCATION") # Repository location

  # Foldernames, you can change your folder names in the quotes
    FAST5="FAST5"           # folder for fast5 raw data
    FASTQ_raw="FASTQ"       # folder for basecalled reads (albacore out)
    FASTQ="DEMULTIPLEXED"   # demultiplex & trimmed output folder (porechop out)
    FASTA_raw="FASTA"       # assembly output folder (assembler out)
    ASSEMBLY="ASSEMBLY"     # auto cp command from assembler out (FASTA_raw) to this
    TAX="TAXONOMY"          # read classification output
    # created at the start
    mkdir -p ${FAST5} ${FASTQ_raw} ${FASTQ} ${FASTA_raw} ${ASSEMBLY} ${TAX}

  # CPU cores
    CPU=$(lscpu -p | egrep -v '^#' | wc -l) # can be changed to e.g. CPU="16"
  # colours, needed for echos
    RED='\033[0;31m'
    GRE='\033[0;32m'
    YEL='\033[0;33m'
    NC='\033[0m'

## PARAMETERS ##
# CONTAINS PARAMETERS FOR PROGRAMS, change at your own risk
  # read length minimum for wtdbg2 assembler
    wtdbg2_readlength="2000"
    wtdbg2_edge_depth="2" # def. 3
    wtdbg2_subsampling="1" # def. 4

###############
### Modules ###
###############
## File download ##

## Basecalling ##
albachore_input()
  {
  echo -e "${RED}Input for Albachore ${NC}"
  echo -e "Enter flowcell type (e.g. ${GRE}FLO-MIN106${NC} or ${GRE}FLO-MIN107${NC}) and hit [Enter]"
  read flowcell
  echo -e "Enter_library kit (e.g. ${GRE}SQK-LSK108${NC}, ${GRE}SQK-RBK004${NC} or ${GRE}SQK-RAD004${NC}) and hit [Enter]"
  read kittype
  }

albachore_execute()
  {
  docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/albacore \
  read_fast5_basecaller.py -r -i /${WORKDIRNAME}/${FAST5} -f $flowcell -t $CPU -q 0 -o fastq -k ${kittype} -r -s ${WORKDIRNAME}/${FASTQ_raw}/
  }

## demultiplexing & trimming ##
porechop_demultiplex()
  {
  docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/porechop \
  porechop -t ${CPU} -i /${WORKDIRNAME}/${FASTQ_raw} -b /${WORKDIRNAME}/${FASTQ}
  # remove untagged barcode
  }

## Assembler ##
wtdbg2_execute()
  {
  echo -e "${RED} Running wtdbg2 with min. read length ${wtdbg2_readlength}bp; -e ${wtdbg2_edge_depth}; -S ${wtdbg2_subsampling}${NC}"
  for fastqfile in ${FASTQ}/*.fastq; do
  filename=$(basename $fastqfile)
    # create assembly map and stuff
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
    wtdbg2 -t $CPU -x ont -e ${wtdbg2_edge_depth} -S ${wtdbg2_subsampling} --rescue-low-cov-edges -i /${WORKDIRNAME}/${fastqfile} -o /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq} # -L ${wtdbg2_readlength}
    # create contigs
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
    wtpoa-cns -t $CPU -i /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq}.ctg.lay.gz -fo /${WORKDIRNAME}/${ASSEMBLY}/${filename%.fastq}.fa
  done
  }
centrifuge_execute()
  {
  for fastqfile in ${FASTQ}/*.fastq; do
    seqID=$(basename ${fastqfile%.fastq})
    docker run --rm -i -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/centrifuge \
    centrifuge -p $CPU -x /centrifuge/database/p_compressed -U /${WORKDIRNAME}/${fastqfile} -S /${WORKDIRNAME}/${TAX}/${seqID}.txt --report-file /${WORKDIRNAME}/${TAX}/${seqID}.log
    # create report for pavian
    docker run --rm -i -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/centrifuge \
    centrifuge-kreport -x /centrifuge/database/p_compressed /${WORKDIRNAME}/${TAX}/${seqID}.txt > ${WORKDIRPATH}/${TAX}/${seqID}_pavian_report.csv
  done
  }

plasflow_execute()
  {
  for fastafile in ${ASSEMBLY}/*.fa; do
    docker run --rm -i -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/plasflow \
    filter_sequences_by_length.pl -input /${WORKDIRNAME}/${fastafile} -output /${WORKDIRNAME}/${fastafile}sta -thresh 2000
    docker run --rm -i -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/plasflow \
    PlasFlow.py --input /${WORKDIRNAME}/${fastafile}sta --output /${WORKDIRNAME}/${TAX}/plasflow_predictions.tsv --threshold 0.7
  done
  }

############################
###   Start of script    ###
############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "
while true; do
    echo -e "${GRE}What do you want to do? [b] [m] [t] [p] or [e]xit${NC}"
    echo -e "[b] ${YEL}Basecalling${NC} - albacore, porechop"
    echo -e "[m] ${YEL}Metagenomassembly${NC} - wtdbg2"
    echo -e "[t] ${YEL}Taxonomy${NC} - centrifuge"
    echo -e "[p] ${YEL}Plasmids${NC} - plasflow"
    read -p "basecalling[b] meta-assembly[m] taxonomy[t] plasmids[p] or [e]xit: " bmtpe
    case $bmtpe in
        [Bb]* ) albachore_input; albachore_execute; porechop_demultiplex; break;;
        [Mm]* ) wtdbg2_execute; break;;
        [Tt]* ) centrifuge_execute; break;;
        [Pp]* ) plasflow_execute; break;;
        [Ee]* ) echo "Exiting script, bye bye"; exit;;
        * ) echo "  Please answer [b] [m] [t] [p] or [e].";;
    esac
done
