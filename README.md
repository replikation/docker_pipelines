# Docker pipelines

![](https://img.shields.io/badge/made%20with-shell-green.svg)
![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)

## Installation
### simple via apt

* install docker via `sudo apt install`
* add the docker group to your user
* execute the appropriate dockers (scripts are inside installations_scripts/)

### better via docker Repository

* either execute the install docker ubuntu script or do it manually
  * it follows the install routine from docker (adding the docker repo to your computer)
* this installs the current docker version for ubuntu, also **adds docker group to user**
* pull all the docker images needed for the scripts
   * scripts are inside installations_scripts/


## Scripts
### cb_docker_ASSEMBLY.sh

### cb_docker_METAGENOME.sh
* if you are using the wtdbg2 assembler adjust the -L flag option in the script to your needs

