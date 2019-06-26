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

  # Serverlocation
    fast5files_server="/volume1/sequencing_data/raw_data" # location of each run_dir with fast5 files
  # IP address USERNAME@IP location in cfg
    IP=$(cat $SCRIPTPATH/cfg) # cfg can also put somewere else (e.g. home)
  # ssh key location
    ssh_key="$HOME/.ssh/id_rsa" # key location
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
basecaller_ask()
  {
  echo -e "${RED}Input for basecaller ${NC}"
  scp ${IP}:${fast5files_server}/${sequence_ID}/flowcell_kit.txt $WORKDIRPATH 2>/dev/null
  cat flowcell_kit.txt 2>/dev/null
  echo -e "Enter flowcell type (e.g. ${GRE}FLO-MIN106${NC} or ${GRE}FLO-MIN107${NC}) and hit [Enter]"
  read flowcell
  echo -e "Enter_library kit (e.g. ${GRE}SQK-LSK108${NC}, ${GRE}SQK-RBK004${NC} or ${GRE}SQK-RAD004${NC}) and hit [Enter]"
  read kittype
  }

albacore_execute()
  {
  docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/albacore \
  read_fast5_basecaller.py -r -i /${WORKDIRNAME}/${FAST5} -f $flowcell -t $CPU \
  -q 0 -o fastq -k ${kittype} -r -s /${WORKDIRNAME}/${FASTQ_raw}/
  }

guppy_cpu()
  {
    mkdir ${FASTQ}
    CPU_half=$(echo $(($CPU / 2)))
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/guppy \
    guppy_basecaller -r -t ${CPU_half} --runners 2 -i /${WORKDIRNAME}/${FAST5} -s /${WORKDIRNAME}/${FASTQ_raw}/ \
    --flowcell ${flowcell} --kit ${kittype} --enable_trimming on --trim_strategy dna -q 0
  }


## demultiplexing & trimming ##
porechop_ask()
  {
  read -p "Are you using barcodes? [yes] or [no] " yn
  case $yn in
      [Yy]* ) barcodes="1"; return;;
      [Nn]* ) barcodes="0"; return;;
      * ) echo "  Please answer [y] or [n].";;
  esac
  }

porechop_execute()
  {
  if (($barcodes==1))
  then
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/porechop \
    porechop -t ${CPU} -i /${WORKDIRNAME}/${FASTQ_raw} -b /${WORKDIRNAME}/${FASTQ}
  else
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/porechop \
    porechop -t ${CPU} -i /${WORKDIRNAME}/${FASTQ_raw} -o /${WORKDIRNAME}/${FASTQ}/all_reads.fastq
  fi
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

## Polishing ##
nanopolish_execute()
  {
  cp ${ASSEMBLY}/*.fa ${POLISH} # copy raw assemblies to polishing sector
  # fastq index; indexed via docker location. should work (only) in a docker mount
  for fastqfile in ${FASTQ}/*.fastq; do
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/nanopolish \
    nanopolish index -s /${WORKDIRNAME}/${FASTQ_raw}/sequencing_summary.txt -d /${WORKDIRNAME}/${FAST5} /${WORKDIRNAME}/${fastqfile}
  done

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
  # polishing
  for assemblyfile in ${POLISH}/*.fa; do
    	sampleID=$(basename ${assemblyfile%.fa})
      mkdir -p ${POLISH}/${sampleID}
    # nanopolish to parallel pipe with variant calling
      docker run --rm -i -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/nanopolish \
      nanopolish_makerange.py /${WORKDIRNAME}/${assemblyfile} |\
      docker run --rm -i -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/nanopolish \
      parallel --results /${WORKDIRNAME}/${POLISH}/${sampleID}/nanopolish.results -P ${CPU} \
      nanopolish variants --consensus -o /${WORKDIRNAME}/${POLISH}/${sampleID}/polished.{1}.vcf -w {1} \
      -r /${WORKDIRNAME}/${FASTQ}/${sampleID}.fastq -b /${WORKDIRNAME}/${assemblyfile%.fa}.sorted.bam -g /${WORKDIRNAME}/${assemblyfile} -t 1 -m 0.1
    # get polished assembly out of results
      docker run --rm -i -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/nanopolish \
      sh -c "nanopolish vcf2fasta -g /${WORKDIRNAME}/${assemblyfile} /${WORKDIRNAME}/${POLISH}/${sampleID}/polished.*.vcf" > ${WORKDIRPATH}/${ASSEMBLY}/${sampleID}_polished.fasta
  done
  }

## other ##
renaming_sequences()
  {
  while IFS=$' ' read barcode SeqID; do
    case "$SeqID" in
    \S.*) mkdir ${ASSEMBLY}/${SeqID}/ ; mv ${ASSEMBLY}/${barcode}_polished.fasta ${ASSEMBLY}/${SeqID}/; cp ${FAST5}/*.txt ${ASSEMBLY}/${SeqID}/ ;;
    *) echo "No SeqID for Barcode $barcode found." ;;
    esac
  done < "${FAST5}/barcodes.txt"
  }

############################
###   Start of script    ###
############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "
while true; do
    echo -e "${GRE}What do you want to do? [u] [p] [a] [n] or [e]${NC}"
    echo -e "[u] ${YEL}UKJ${NC} - fast5 download, guppy_cpu, porechop, wtdbg2, nanopolish, rename"
    echo -e "[p] ${YEL}Pipeline${NC} - guppy_cpu, porechop, wtdbg2, nanopolish"
    echo -e "[a] ${YEL}ASSEMBLY${NC} - porechop, wtdbg2, fastq files in ${FASTQ_raw}/"
    echo -e "[n] ${YEL}POLISH${NC} - nanopolish, .fast5 in ${FAST5}/, .fastq in ${FASTQ}/, .fa in ${ASSEMBLY}"
    read -p "UKJ[f] Pipeline[p] assembly[a] polish[n] exit[e]: " upante
    case $upante in
        [Uu]* ) Downloadinput; basecaller_ask; porechop_ask; Download_FAST5; guppy_cpu; porechop_execute; wtdbg2_execute; nanopolish_execute; renaming_sequences; break;;
        [Pp]* ) basecaller_ask; porechop_ask; guppy_cpu; porechop_execute; wtdbg2_execute; nanopolish_execute; break;;
        [Aa]* ) porechop_ask; porechop_execute; wtdbg2_execute; break;;
        [Nn]* ) nanopolish_execute; break;;
        [Tt]* ) renaming_sequences; break;;
        [Ee]* ) echo "  Exiting script, bye bye"; exit;;
        * ) echo "  Please answer [f] [a] [n] or [e].";;
    esac
done
