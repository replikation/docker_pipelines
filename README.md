# Docker pipelines

![](https://img.shields.io/badge/made%20with-shell-green.svg)
![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)

## Installation
### Docker

* install `docker` via `sudo apt install docker.io` 
	* or `docker-ce` via [this](https://docs.docker.com/install/linux/docker-ce/ubuntu/)
* add the docker group to your user via `sudo usermod -a -G docker $USER`
* reboot after this

### Script-options
* Scripts are getting all neccessary components automatically if missing
* all scripts have "options", open them via editor and check the option section

## Scripts
### cb_docker_ASSEMBLY.sh

### cb_docker_METAGENOME.sh
* if you are using the wtdbg2 assembler adjust the -L flag option in the script to your needs

### cb_docker_ILLUMINA.sh
