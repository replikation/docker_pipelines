#!/usr/bin/env bash
echo "Stoping and removing all dangling images"
docker stop $(docker ps -a -q)
docker rm $(docker ps -a -q)
docker rmi $(docker images -f "dangling=true" -q)
