#!/usr/bin/env bash

echo "pulling docker images"
docker pull replikation/guppy
docker pull replikation/porechop
docker pull replikation/wtdbg2
docker pull replikation/centrifuge
docker pull replikation/plasflow
echo "removing dangling images"
docker rmi $(docker images -f "dangling=true" -q)
