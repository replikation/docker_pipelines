#!/usr/bin/env bash
SCRIPT=$(readlink -f "$0") #  e.g. /home/user/bin/foo.sh
SCRIPTPATH=$(dirname "$SCRIPT") # e.g. /home/user/bin

# outputfolder for everything albachore related
mkdir -p albachore_output
## Options
# flowcell type
  echo -e "Enter flowcell type (e.g. ${GRE}FLO-MIN106${NC} or ${GRE}FLO-MIN107${NC}) and hit [Enter]"
  read flowcell
    echo "$flowcell" > ${SCRIPTPATH}/flowcell_used.txt
# kittype
  echo -e "Enter_library kit (e.g. ${GRE}SQK-LSK108${NC}, ${GRE}SQK-RBK004${NC} or ${GRE}SQK-RAD004${NC}) and hit [Enter]"
  read kittype
    echo "$kittype" > ${SCRIPTPATH}/kittype_used.txt
