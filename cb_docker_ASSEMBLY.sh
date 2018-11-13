#!/usr/bin/env bash
#!/bin/bash
#!/usr/bin/bash

# to do
# check the output of porechop for unbarcoded sequences (are they also called BC or none?)

## Docker ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found. Please run the installation skript, Aborting.${NC}"; exit 1; }
  echo "Docker identified"

## DIR LOCATIONS ##
  WORKDIRPATH=$(pwd)
  WORKDIRNAME=${PWD##*/}
  SCRIPTLOCATION=$(readlink -f "$0")
  SCRIPTPATH=$(dirname "$SCRIPTLOCATION")

  # Foldernames
  FAST5="FAST5"         # folder for fast5 raw data
  FASTQ_raw="FASTQ"     # folder for basecalled fastq (alba out)
  FASTQ="DEMULTIPLEXED" # demultiplex output folder (porech. out)
  FASTA_raw="FASTA"     # assembly output folder (assembler out)
  ASSEMBLY="ASSEMBLY"   # standard location for assembly (standard output)
  POLISH="POLISHING"    # polish output folder

  mkdir -p ${FAST5} ${FASTQ_raw} ${FASTQ} ${FASTA_raw} ${POLISH} ${ASSEMBLY}

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
  YEL='\033[0;33m'
  NC='\033[0m'

## PARAMETERS
  # read length minimum for wtdbg2 assembler
  wtdbg2_readlength="5000"

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
  scp ${IP}:${fast5files_server}/${sequence_ID}/flowcell_kit.txt $WORKDIRPATH 2>/dev/null
  cat flowcell_kit.txt 2>/dev/null
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
  rm ${FASTQ}/none.fastq 2>/dev/null
  }

## Assembler ##
wtdbg2_execute()
  {
  for fastqfile in $FASTQ/*.fastq; do
  filename=$(basename $fastqfile)
    #create assembly map and stuff
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
    wtdbg2 -t $CPU -i /${WORKDIRNAME}/${fastqfile} -o /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq} -L ${wtdbg2_readlength}
    #create contigs
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
    wtpoa-cns -t $CPU -i /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq}.ctg.lay.gz -fo /${WORKDIRNAME}/${ASSEMBLY}/${filename%.fastq}.fa
  done
  }

## Polishing ##
nanopolish_execute()
{
  cp ${ASSEMBLY}/*.fa ${POLISH} # copy raw assemblies to polishing sector
##UNTESTED
  # fastq index
  #for fastqfile in ${FASTQ}/*.fastq; do
  #  docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/nanopolish \
  #  nanopolish index -s /${WORKDIRNAME}/${FASTQ}/sequencing_summary.txt -d /${WORKDIRNAME}/${FAST5} /${WORKDIRNAME}/${fastqfile}
  #done

  # align against metagenome
  for assemblyfile in ${POLISH}/*.fa ; do
      filename=$(basename $assemblyfile)
      # bwa index
        docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/bwa \
    	  bwa index /${WORKDIRNAME}/$assemblyfile
      # bwa align
        docker run --rm -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/bwa \
        bwa mem -x ont2d -t $CPU /${WORKDIRNAME}/$assemblyfile /${WORKDIRNAME}/${FASTQ}/${filename%.fa}.fastq > /${WORKDIRPATH}/${assemblyfile%.fa}.sam
      # sam conversion
        docker run --rm -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/samtools \
        samtools view -Sb /${WORKDIRNAME}/${assemblyfile%.fa}.sam  > /${WORKDIRPATH}/${assemblyfile%.fa}.bam
      # bam sort
        docker run --rm -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/samtools \
        samtools sort -@$CPU /${WORKDIRNAME}/${assemblyfile%.fa}.bam -T /tmp/aln.sorted -o /${WORKDIRNAME}/${assemblyfile%.fa}.sorted.bam
      # bam index
        docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/samtools \
        samtools index /${WORKDIRNAME}/${assemblyfile%.fa}.sorted.bam
   done
##UNTESTED
  # polishing
    ##python nanopolish_makerange.py draft.fa | parallel --results nanopolish.results -P 8 \
    ##nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r reads.fa -b reads.sorted.bam -g draft.fa -t 4 --min-candidate-frequency 0.1
  # put together
    ## nanopolish vcf2fasta -g draft.fa polished.*.vcf > ${WORKDIRPATH}/${POLISH}/polished_genome.fasta
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
    echo -e "[a] ASSEMBLY - put one fastq file for each sample in ${YEL}DEMULTIPLEX${NC}"
    echo "[t] for testing modules (ignore for normal usage)"
    read -p "FULL[f] metagenome[m] assembly_only[a] nanopolish[n] exit[e]: " fmante
    case $fmante in
        [Ff]* ) Downloadinput; albachore_input; Download_FAST5; albachore_execute; porechop_demultiplex; wtdbg2_execute; break;;
        [Mm]* ) placeholder; break;;
        [Aa]* ) wtdbg2_execute; break;;
        [Nn]* ) nanopolish_execute; break;;
        [Tt]* ) albachore_input; albachore_execute; break;;
        [Ee]* ) echo "  Exiting script, bye bye"; exit;;
        * ) echo "  Please answer [f] [m] [a] [n] or [e].";;
    esac
done
