#!/usr/bin/env bash
#!/bin/bash
#!/usr/bin/bash

## Docker ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found. Please run the installation skript, Aborting.${NC}"; exit 1; }
  echo "Docker identified"

## DIR LOCATIONS ##
  WORKDIRPATH=$(pwd)
  WORKDIRNAME=${PWD##*/}
  SCRIPTLOCATION=$(readlink -f "$0")
  SCRIPTPATH=$(dirname "$SCRIPTLOCATION")

  # Foldernames
  FAST5="FAST5_raw"
  FASTQ_raw="FASTQ_albachore"
  FASTQ="DEMULTIPLEXED_porechop"
  FASTA_raw="ASSEMBLY_wtdbg2"
  FASTA="POLISHING_nanopolish"
  mkdir -p $FAST5 $FASTQ_raw $FASTQ $FASTA_raw $FASTA

  # Serverlocation
  fast5files_server="/volume1/sequencing_data/raw_data"
  # IP address USERNAME@IP location in cfg
  IP=$(cat $SCRIPTPATH/cfg)
  ssh_key="$HOME/.ssh/id_rsa"
  # CPU cores
  CPU=$(lscpu -p | egrep -v '^#' | wc -l)

  # colours
  RED='\033[0;31m'
  #BLU='\033[0;34m'
  GRE='\033[0;32m'
  #YEL='\033[0;33m'
  NC='\033[0m'

## PARAMETERS
  # read length minimum for wtdbg2 assembler
  wtdbg2_readlength="7000"

###############
### Modules ###
###############
## File download ##
Downloadinput()
  {
  echo "The following Sequence runs were found:" ; echo " "
  ssh $IP ls $fast5files_server
  echo " "
  echo -e "Please enter foldername for download (e.g. ${GRE}2018.04.19.iimk_run${NC}) and hit [ENTER]"
  read sequence_ID
  scp ${IP}:${fast5files_server}/${sequence_ID}/flowcell_kit.txt $WORKDIRPATH
  }

Download_FAST5()
  {
  read -p "Do you want to download this folder: $sequence_ID  [yes/no]: " yn
  case $yn in
      [Yy]* ) echo -e "${RED}Starting download ${NC}"; rsync --rsync-path=/bin/rsync  -r -e "ssh -i $ssh_key" $IP:${fast5files_server}/${sequence_ID}/* $FAST5 ; echo -e "${RED}Download done ${NC}";;
      [Nn]* ) Download_FAST5 ;;
      * ) echo "Please answer yes or no.";;
  esac
  }

## Basecalling ##
albachore_input()
  {
  echo -e "${RED}Input for Albachore ${NC}"
  cat flowcell_kit.txt 2>/dev/null
  echo -e "Enter flowcell type (e.g. ${GRE}FLO-MIN106${NC} or ${GRE}FLO-MIN107${NC}) and hit [Enter]"
  read flowcell
  echo -e "Enter_library kit (e.g. ${GRE}SQK-LSK108${NC}, ${GRE}SQK-RBK004${NC} or ${GRE}SQK-RAD004${NC}) and hit [Enter]"
  read kittype
  }

albachore_execute()
  {
  docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/albacore \
  read_fast5_basecaller.py -r -i $WORKDIRPATH/$FAST5 -f $flowcell -t $CPU -q 0 -o fastq -k $kittype -r -s $WORKDIRPATH/$FASTQ_raw/
  }

## demultiplexing & trimming ##
  # not tested, for unbarcoded samples
porechop_demultiplex()
  {
  docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/porechop \
  porechop -t $CPU -i /${WORKDIRNAME}/${FASTQ_raw}/ -b /${WORKDIRNAME}/${FASTQ}
  # remove untagged barcode
  rm ${FASTQ}/none.fastq 2>/dev/null
  }

## Assembler ##
wtdbg2_execute()
  {
  for fastqfile in $FASTQ/*.fastq; do
  filename=$(basename $fastqfile)
  echo "$filename"
    echo "$filename assembly"
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
    wtdbg2 -t $CPU -i /${WORKDIRNAME}/${FASTQ}/${filename} -o /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq} -L $wtdbg2_readlength
  done
  }

## Polishing ##
nanopolish_execute()
{
  # fastq index
  for x in $FASTQ/*.fastq; do
    fastqfile=$(basename $x)
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/nanopolish \
    nanopolish index -d ${WORKDIRNAME}/$FAST5 ${WORKDIRNAME}/$FASTQ/$fastqfile
  done


  # align against metagenome
    # for x in assembly draft

    #bwa index draft.fa #location of assemblys unclear

    #bwa mem -x ont2d -t 8 draft.fa reads.fa | samtools sort -o reads.sorted.bam -T reads.tmp -
    #samtools index reads.sorted.bam
  # polishing
    #python nanopolish_makerange.py draft.fa | parallel --results nanopolish.results -P 8 \
      #nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r reads.fa -b reads.sorted.bam -g draft.fa -t 4 --min-candidate-frequency 0.1
  # put together
    #nanopolish vcf2fasta -g draft.fa polished.*.vcf > polished_genome.fa
}

placeholder()
{
  echo "placeholder"
}

############################
### Start of script      ###
############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "
while true; do
    echo -e "${GRE}What do you want to do? [f] [m] [a] [n] or [e]${NC}"
    echo "[f] FULL - fast5 download, basecalling, demultiplex, assembly, polishing"
    echo "[t] for testing modules (ignore for normal usage)"
    read -p "FULL[f] metagenome[m] assembly_only[a] nanopolish[n] exit[e]: " fmante
    case $fmante in
        [Ff]* ) Downloadinput; albachore_input; Download_FAST5; albachore_execute; porechop_demultiplex; wtdbg2_execute; break;;
        [Mm]* ) placeholder; break;;
        [Aa]* ) wtdbg2_execute; break;;
        [Nn]* ) placeholder; break;;
        [Tt]* ) Downloadinput; break;;
        [Ee]* ) echo "  Exiting script, bye bye"; exit;;
        * ) echo "  Please answer [f] [m] [a] [n] or [e].";;
    esac
done
