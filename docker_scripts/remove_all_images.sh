#!/usr/bin/env bash
echo "Removing all Docker images now"

docker rmi -f $(docker images -a -q)
