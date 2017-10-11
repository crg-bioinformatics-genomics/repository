#!/bin/bash
# bash RBProfiles.sh <random_number> <email> <title> <proteinFile>
echo $1 $2 $3 $4

set -e
#set -o pipeline

cp -r template tmp/$1
cd   tmp/$1

	cp "$4" data/inseq.fasta

	#python RBPprofilerBeta.py 2&>log
	python RBPprofilerGamma.py 2&>log

cd ..
