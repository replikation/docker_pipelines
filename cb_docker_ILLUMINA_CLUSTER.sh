#!/usr/bin/env bash

## Docker ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found. Please run the installation skript, Aborting.${NC}"; exit 1; }

## OPTIONS ##
  WORKDIRPATH=$(pwd) # for docker mountpoint (-v)
  WORKDIRNAME=${PWD##*/} # for docker mountpoint (-v)

  # FASTQ prefixes for forward and reverse pairs, change the quotes
    FPREFIX="_R1_" # e.g. FPREFIX="_1.fasta"
    RPREFIX="_R2_" # e.g. FPREFIX="_2.fasta"

  # Foldernames, you can change your folder names in the quotes
    # Server locations (for option a - the complete automated workflow)
    FASTQ_serv="/volume1/Shared_Regensburg/FASTQ_Rohdaten"        # RAWDATA on Server
    ASSEMBLY_serv="/volume1/Shared_Regensburg/ASSEMBLIES_fertig"  # ASSEMBLY on Server
    CLUSTER_serv="/volume1/Shared_Regensburg/ANALYSEN"            # CLUSTER Results on Server
    IP=$(cat ~/cfg)                                               # location of USERNAME@IPADRESS
      # make sure you work with ssh keys!
    # Foldernames in Current Workdir (for all options)
    FASTQ_Wdir="FASTQ"                  # FASTQ location (serv & Wdir)
    ASSEMBLY="ASSEMBLY"                 # Results of Assemblies
    CLUSTER="CLUSTER"                   # Results cluster analysis

  # CPU cores
    CPU=$(lscpu -p | egrep -v '^#' | wc -l) # can be changed to e.g. CPU="16"
  # colours, needed for echos
    RED='\033[0;31m'
    GRE='\033[0;32m'
    YEL='\033[0;33m'
    NC='\033[0m'

###############
### Modules ###
###############
serv_check_downl()
  {
    echo "Download & Assembly files from $FASTQ_serv ?"
    read -p "Proceed? [yes] or [no] " yn
    case $yn in
        [Yy]* ) echo "Starting Pipeline";;
        [Nn]* ) echo "Exiting..."; exit;;
        * ) echo "  Please answer [y] or [n].";;
    esac
    RAW_files=$(ssh $IP ls ${FASTQ_serv})
    while IFS= read -r samplename ; do
       echo "$samplename"
       if ssh < /dev/null $IP stat ${ASSEMBLY_serv}/${samplename} \> /dev/null 2\>\&1; then
            echo "${samplename} exists"
          else
            echo "${samplename} does not exist, starting download"
            mkdir -p $FASTQ_Wdir
            scp -r ${IP}:${FASTQ_serv}/${samplename} ${FASTQ_Wdir}
       fi
    done < <(printf '%s\n' "$RAW_files")
  }

## Read QC ##
read_quality()
  {
    mkdir -p $ASSEMBLY
    for illuminarun in ${FASTQ_Wdir}/* ; do
      F_file=$(ls ${illuminarun}/*${FPREFIX}*)
      R_file=$(ls ${illuminarun}/*${RPREFIX}*)
      SampleID=${illuminarun##*/}
      mkdir -p ${ASSEMBLY}/${SampleID}
      # fastqc run
      docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/fastqc \
      ${WORKDIRNAME}/${F_file} /${WORKDIRNAME}/${R_file} -t $CPU --outdir=/${WORKDIRNAME}/${ASSEMBLY}/${SampleID}
    done
    rm -f ${ASSEMBLY}/*/*.zip
  }

## Assembler ##
unicycler_execute()
  {
    for illuminarun in ${FASTQ_Wdir}/* ; do
      F_file=$(ls ${illuminarun}/*${FPREFIX}*)
      R_file=$(ls ${illuminarun}/*${RPREFIX}*)
      SampleID=${illuminarun##*/}
      mkdir -p ${ASSEMBLY}/${SampleID}
      # fastqc run
       docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/unicycler \
       -1 /${WORKDIRNAME}/$F_file -2 /${WORKDIRNAME}/$R_file -t $CPU -o /${WORKDIRNAME}/${ASSEMBLY}/${SampleID}
      # rename assembly file
      mv ${ASSEMBLY}/${SampleID}/assembly.fasta ${ASSEMBLY}/${SampleID}/${SampleID}.fasta
    done
    # remove all gfa files
    rm -f ${ASSEMBLY}/*/*.gfa
  }

## sourmash ##
serv_get_fasta()
  {
    # this is to supplement the cluster data with more fasta's
    echo "Downloading Assemblies from Server to temporary signature folder"
    mkdir -p tmp_signatures
    scp -r ${IP}:${ASSEMBLY_serv}/*/*.fasta tmp_signatures
  }

sourmash_execute()
  {
  # creating signatures
  echo "Starting sourmash, it breaks if two or more seq.IDs are identical."
  mkdir -p tmp_signatures
  cp ${ASSEMBLY}/*/*.fasta tmp_signatures/
  for fasta in tmp_signatures/*.fasta; do mv ${fasta} ${fasta%.fasta}; done
  cd tmp_signatures/
  for sequence in *; do
    docker run --rm -it -v ${WORKDIRPATH}/tmp_signatures:/tmp_signatures replikation/sourmash \
    /bin/sh -c "cd tmp_signatures/ && sourmash compute -n 5000 -k 31,51 ${sequence} -o ${sequence}.sig"
  done
  cd ..
  # create index out of signatures
  docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/sourmash \
  /bin/sh -c "sourmash compare /${WORKDIRNAME}/tmp_signatures/*.sig -o /${WORKDIRNAME}/tmp_signatures/index_sig"
  echo "comparing done"
  # plotting files to pdf
  docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/sourmash \
  /bin/sh -c "sourmash plot --pdf --labels /${WORKDIRNAME}/tmp_signatures/index_sig && mv *.pdf /${WORKDIRNAME}/${CLUSTER}"
  # remove temporary signature folder
  rm -rf tmp_signatures
  }

## file upload ##
serv_upload()
  {
    echo "Uploading Cluster and Assembly files to server"
    scp -r ${ASSEMBLY}/* $IP:$ASSEMBLY_serv #check how the file structure looks after upload
    scp -r ${CLUSTER}/* $IP:$CLUSTER_serv #check how the file structure looks after upload
  }

############################
###   Start of script    ###
############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "
while true; do
    echo "Illumina Paired end, automated Pipeline"
    echo -e "The following fastq IDs for forward/reverse are used: ${YEL}${FPREFIX}${NC} and ${YEL}${RPREFIX}${NC}"
    echo " "
    echo -e "${GRE}What do you want to do? [a] [c] [s] [e]xit${NC}"
    echo -e "[a] ${YEL}Automated server analysis - assembly, clusteranalysis${NC}"
    echo "    Start in a empty folder, downloads and uploades files automatically"
    echo -e "[c] ${YEL}Current Workingdir Assembly${NC} - unicycler assembly"
    echo "    Expects one folder for each illumina run in ./$FASTQ_Wdir/, e.g $FASTQ_Wdir/e.coli_1 $FASTQ_Wdir/k.pneu"
    echo -e "[s] ${YEL}Current Workingdir Clusteranalysis${NC} - sourmash cluster analysis"
    echo "    Expects one folder + .fasta for each assembly in ./$ASSEMBLY/, e.g $ASSEMBLY/e.coli_1/e.coli_1.fasta "
    echo " "
    read -p "Automated[a] Current Wdir Assembly[c] Sourmash cluster [s] or [e]xit: " acsde
    case $acsde in
        [Aa]* ) serv_check_downl; read_quality; unicycler_execute; serv_get_fasta; sourmash_execute; serv_upload; break;;
        [Cc]* ) read_quality; unicycler_execute; break;;
        [Ss]* ) sourmash_execute; break;;
        [Dd]* ) serv_upload; break;;
        [Ee]* ) echo "Exiting script, bye bye"; exit;;
        * ) echo "  Please answer [a] [c] [s] or [e].";;
    esac
done
echo -e "${GRE} Finished all steps ${NC}"
