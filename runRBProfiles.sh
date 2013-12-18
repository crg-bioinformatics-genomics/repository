#!/bin/bash
# bash RBProfiles.sh <random_number> <title> <proteinFile>
echo $1 $2 $3

set -e
#set -o pipeline

cp -r template tmp/$1
cd   tmp/$1  

	python ...

	cp ./matrix.png  ./outputs/Heatmap.$1.$3.png

cd ..
