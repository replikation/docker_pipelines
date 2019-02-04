#!/usr/bin/env bash
echo "removing dangling images"
docker stop $(docker ps -a -q)
docker rmi $(docker images -f "dangling=true" -q)
