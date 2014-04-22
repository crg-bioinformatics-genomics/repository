from __future__ import division
import os
import sys
import random
import subprocess
import pickle
import math
import numpy

#sort -nk3 ROCminimumError.txt | awk '($4!="is"){print $1" "$2" "$4}' > Error_motif_scale_threshold.txt
#sort -nk3 ROCmaxCorr.txt | uniq | awk '($4!="threshold"&& $4>0.6){print $1" "$2" "$4}' > Corr_Domain_scale_threshold.txt

#ENSP00000428756_9-72 66 0.860546875 domain+scale
#subprocess.call("cat /home/carmen/Documents/CRG/Profile_motifs/new2014/results/DomainsAndBestScale_pval0.1.txt",shell=True)
word= "awk \'{ print $1\" \"$2}\' /home/carmen/Documents/CRG/Profile_motifs/new2014/results/DomainsAndBestScale_pval0.1.NEW > domain_scale.txt"
#print word
#let=raw_input(": ")
subprocess.call(word,shell=True)
#grep -f domain_scale.txt ROCminimumError.txt > j1.txt"
#grep -f domain_scale.txt ROCaccumAndThresSubset0.txt  

fileA=open("domain_scale.txt","r")
dosca={}
for line in fileA:
	dosca[line[:-1]]=1
print dosca
#ROCminimumError.txt
#domain scale auc bestThres sens spec

fileA=open("ROCminimumErrorSubset0.txt","r")
final={}
for line in fileA:
	if line[0]!="#" and line.split()[3]!="is" and line.split()[4]!="threshold":
		print line
		coppia=line.split()[0]+" "+line.split()[1]
		if coppia in dosca:
			final[coppia]=line[:-1]
print final



#corr domain scales auc threshold spec sens

fileA=open("ROCaccumAndThresSubset0.txt","r")
bestCorr={}
for line in fileA:
	if line[0]!="#" and line.split()[3]!="is" and line.split()[4]!="threshold":
		coppia=line.split()[1]+" "+line.split()[2]
		if coppia in dosca:
			if coppia in bestCorr:
				if float(bestCorr[coppia].split()[3])<=float(line.split()[3]):
					bestCorr[coppia]=line
			else:
				bestCorr[coppia]=line
outfile=open("joinedCorrError.NEW","w")
outfile.write("#domain scales aucError bestThresError aucCorr bestCorr thresholdForAccum\n")
for key in final:
	outfile.write(key+" "+str(round(float(final[key].split()[2]),2))+" "+final[key].split()[3]+" "+str(round(float(bestCorr[key].split()[3]),2))+" "+bestCorr[key].split()[0]+" "+bestCorr[key].split()[4]+"\n")

