#!/usr/bin/env bash
#!/bin/bash
#!/usr/bin/bash

### README ###
# This will be now the new pipeline using docker containers as programs
# I need to create controlled output dirs for every step
# each docker mounts the current directory when i start this script so it can use the files

### Path's ##
  #currDir=$(pwd)
  SCRIPT=$(readlink -f "$0") #  e.g. /home/user/bin/foo.sh
  SCRIPTPATH=$(dirname "$SCRIPT") # e.g. /home/user/bin
### Parameters & colours ##
  #IP=$(cat $SCRIPTPATH/cfg)
  #CPU=$(lscpu -p | egrep -v '^#' | wc -l)
    RED='\033[0;31m'
    #BLU='\033[0;34m'
    #GRE='\033[0;32m'
    #YEL='\033[0;33m'
    NC='\033[0m' # No Color
### Docker start up ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}  docker not found. Aborting.${NC}"; exit 1; }
  echo "  docker identified"
  #docker pull #### BASIC-debian-DOCKER ##########

### Modules ###
albachore_input()
  {
  bash ${SCRIPTPATH}modules/albachore_input.sh
  }

albachore_execute()
  {
  bash ${SCRIPTPATH}modules/albachore_run.sh
  }

## Assembler ##

wtdbg_execute()
{
 echo "docker run"
}

placeholder()
{
  echo "placeholder"
}


## How it should work
# i only need a docker dependency check
# each program should be called via a module with unified Output
# check hadrien to see how his docker installs look like


############################
# Start of script OLD#
############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "
while true; do
    echo "What do you want to do? [f] [m] [a] [n] or [e]"
    echo "dev: [t] for testing modules"
    read -p "full_pipeline[f] metagenome[m] assembly_only[a] nanopolish[n] exit[e]: " fmante
    case $fmante in
        [Ff]* ) placeholder; break;;
        [Mm]* ) placeholder; break;;
        [Aa]* ) placeholder; break;;
        [Nn]* ) placeholder; break;;
        [Tt]* ) placeholder; break;;
        [Ee]* ) echo "  Exiting script, bye bye"; exit;;
        * ) echo "  Please answer [f] [m] [a] [n] or [e].";;
    esac
done
