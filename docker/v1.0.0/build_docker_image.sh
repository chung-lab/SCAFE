#!/bin/bash
cwd=`dirname $0`
cd $cwd

docker build -t cchon/scafe:v1.0.0 .
docker tag cchon/scafe:v1.0.0 cchon/scafe:latest