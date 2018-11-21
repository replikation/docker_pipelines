#!/usr/bin/env bash

## Docker ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found. Please run the installation skript, Aborting.${NC}"; exit 1; }
  echo "Docker identified"

## OPTIONS ##
  WORKDIRPATH=$(pwd) # for docker mountpoint (-v)
  WORKDIRNAME=${PWD##*/} # for docker mountpoint (-v)
  SCRIPTLOCATION=$(readlink -f "$0") # Script location
  SCRIPTPATH=$(dirname "$SCRIPTLOCATION") # Repository location

  # Foldernames, you can change your folder names in the quotes
    FAST5="FAST5"           # folder for fast5 raw data
    FASTQ_raw="FASTQ"       # folder for basecalled reads (albacore out)
    FASTQ="DEMULTIPLEXED"   # demultiplex & trimmed output folder (porechop out)
    FASTA_raw="FASTA"       # assembly output folder (assembler out)
    ASSEMBLY="ASSEMBLY"     # Assembly (add a cp command to a assembler to transfer the fasta here)
    POLISH="POLISHING"      # polish output folder, copies the polished fasta also to ASSEMBLY
    # created at the start
    mkdir -p ${FAST5} ${FASTQ_raw} ${FASTQ} ${FASTA_raw} ${POLISH} ${ASSEMBLY}

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
    wtdbg2_readlength="8000"

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
  for fastqfile in ${FASTQ}/*.fastq; do
  filename=$(basename $fastqfile)
    #create assembly map and stuff
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
    wtdbg2 -t $CPU -i /${WORKDIRNAME}/${fastqfile} -o /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq} -L ${wtdbg2_readlength}
    #create contigs
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
    wtpoa-cns -t $CPU -i /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq}.ctg.lay.gz -fo /${WORKDIRNAME}/${ASSEMBLY}/${filename%.fastq}.fa
  done
  # renaming the "none assembly" to avoid polishing it (usually bwa gets stuck)
    mv ${FASTQ}/none.fa ${FASTQ}/none.fasta 2>/dev/null
  # remove empty assembly files
    find ${WORKDIRPATH}/${ASSEMBLY}/ -size  0 -print0 |xargs -0 rm --
  }
centrifuge_execute()
  {
  echo "test"
  }

plasflow_execute()
  {
  echo "test"
  }


############################
###   Start of script    ###
############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "
while true; do
    echo -e "${GRE}What do you want to do? [b] [t] [p] or [e]xit${NC}"
    echo -e "[b] ${YEL}Basecalling${NC} - albacore, porechop"
    echo -e "[t] ${YEL}Taxonomy${NC} - centrifuge"
    echo -e "[p] ${YEL}Plasmids${NC} - plasflow"
    read -p "UKJ[f] Pipeline[p] assembly[a] polish[n] exit[e]: " ptb
    case $ptb in
        [Pp]* ) albachore_input; albachore_execute; porechop_demultiplex; break;;
        [Tt]* ) centrifuge_execute; break;;
        [Bb]* ) plasflow_execute; break;;
        [Ee]* ) echo "  Exiting script, bye bye"; exit;;
        * ) echo "  Please answer [f] [a] [n] or [e].";;
    esac
done
