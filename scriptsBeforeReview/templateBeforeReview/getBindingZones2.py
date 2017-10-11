from __future__ import division
import os
import sys
import random
import subprocess
import pickle
import math
import numpy
import scipy

#takes the highest peaks does NOT force the extension and forward only the peaks longer than 50amino acids

infile=open("./data/inseq.fasta","r")
sequence = infile.readlines()[1]


#infile=open("./outputs/profile.txt","r")
infile=open("./outputs/smoothy21.txt","r")
tot=len(infile.readlines())

zones = []

i=0
infile=open("./outputs/smoothy21.txt","r")
for line in infile:
	print line
	#zones[0,i]=round(float(line.split()[2]),2)
	zones.append(round(float(line.split()[1]),2))
	i+=1
infile.close()



maxi=max(zones)
centro=zones.index(maxi)
o=open("./outputs/list.txt","w")

#take peaks 
while (maxi>=0.5):
	maxi=max(zones)
	centro=zones.index(maxi)

	extend=1
	i=0
	piece=[]
	piece.append(maxi)
	stretch_start=centro
	stretch_end=centro
	while(extend==1):
		i+=1
		if (centro-i)>0:
			if zones[centro-i]>=0.4: #if pos zone than extend the stretch
				extend=1
				piece.insert(0,zones[centro-i])
				zones[centro-i]=-100
				stretch_start=centro-i
			else:
				extend=0
		#print len(piece)

	extend=1
	i=0
	while(extend==1):
		i+=1
		if (centro+i)<len(zones):
			if zones[centro+i]>=0.4:
				extend=1
				piece.append(zones[centro+i])
				zones[centro+i]=-100
				stretch_end=centro+i
			else:
				extend=0
		#print len(piece)
	zones[centro]=-100

	print len(piece),"stretch_start ",stretch_start, " stretch_end ",stretch_end
	if (stretch_end-stretch_start)>=50:
		substring=sequence[stretch_start:stretch_end]
		o.write(">zone_"+str(stretch_start)+"_"+str(stretch_end)+"\n")
		o.write(substring+"\n")
	maxi=max(zones)
	centro=zones.index(maxi)

o.close()


