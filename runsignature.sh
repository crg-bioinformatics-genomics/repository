#!/bin/bash
# bash RBProfiles.sh <random_number> <email> <title> <proteinFile>
echo $1 $2 $3 $4

set -e
#set -o pipeline

cp -r template tmp/$1
cd   tmp/$1

	cp "$4" data/inseq.fasta
	NUM_LINES=$(grep ">" data/inseq.fasta | wc -l | awk '{print $1}')
	echo NUM_LINES $NUM_LINES

	if [ $NUM_LINES = 1 ]; then
		echo execution of Gamma > sh.log
		#python RBPprofilerBeta.py 2&>log
		python RBPprofilerGamma.py 2&>log
	else
		echo execution of Multi >sh.log
		#python RBPprofilerMulti.py 2&>log
		python RBPprofilerMulti.py $1 $2 2&>log
	fi


	

cd ..
