#!/usr/bin/env bash
echo "removing dangling images"
docker rmi $(docker images -f "dangling=true" -q)
