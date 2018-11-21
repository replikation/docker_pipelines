#!/usr/bin/env bash

echo "pulling docker images"
docker pull replikation/albacore
docker pull replikation/porechop
docker pull replikation/wtdbg2
docker pull replikation/bwa
docker pull replikation/samtools
docker pull replikation/nanopolish
echo "removing dangling images"
docker rmi $(docker images -f "dangling=true" -q)
