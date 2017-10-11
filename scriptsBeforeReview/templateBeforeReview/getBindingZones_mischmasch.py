from __future__ import division
import os
import sys
import random
import subprocess
import pickle
import math
import numpy
import scipy

####takes all the peaks >=0.5 and forces the extension to at least 50amino acids

infile=open("./data/inseq.fasta","r")
sequence = infile.readlines()[1]

seq=[]
for i in sequence:
	seq.append(i)
seq.remove((seq[len(seq)-1]))


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

o=open("./outputs/list.txt","w")

ll=[]
infile=open("./outputs/HMMER.results","r")
for line in infile:
	ll.append(line)
	

for element in ll:
	print element
	if (int(element.split()[1])-int(element.split()[2])<=54):#### for FMRP iso7  meinetwegen bitte weg damit
		for i in range(int(element.split()[1])-6,int(element.split()[2])+1):
			seq[i]="-"
	else:
		for i in range(int(element.split()[1]),int(element.split()[2])+1):
			seq[i]="-"

HMMERlist=[]
start=0
for i in range(len(seq)):
	if seq[i]=="-" and start==0:
		start=i
	if seq[i]!="-" and start!=0:
		diff=(i-1-start)
		if diff>=50:
			print "OK"
			HMMERlist.append(str(start)+" "+str(i-1))
			o.write(">HMMERzone_"+str(start)+"_"+str(i-1)+"\n")
			o.write(sequence[start:i-1]+"\n")				
		else:
			print "NOTOK"
			extension=50-diff
			print extension 
			newstart=start-int(extension/2)
			newend=i-1+int(extension/2)
			HMMERlist.append(str(newstart)+" "+str(newend))			
			o.write(">HMMERzone_"+str(newstart)+"_"+str(newend)+"\n")
			o.write(sequence[newstart:newend]+"\n")				

		start=0


otherZones=[]
legalZones=[]
start=0
for i in range(len(seq)):
	if seq[i]!="-" and start==0:
		start=i
	if seq[i]=="-" and start!=0:
		subscore=zones[start:i-1]
		otherZones.append(str(start)+" "+str(i-1))			
		print str(start)+" "+str(i-1)
		print (max(subscore))
		if max(subscore)>=0.55:
			diff=i-start
			print i,start,"diff",diff
			if diff<50:
				extension=50-diff
				newstart=start-int(extension/2)
				newend=i-1+int(extension/2)
				legalZones.append(str(newstart)+" "+str(newend))			
			else:
				legalZones.append(str(start)+" "+str(i-1))			

		start=0

for ii in legalZones:
	print ii
	stretch_start=int(ii.split()[0])
	stretch_end=int(ii.split()[1])
	print (stretch_end-stretch_start),"stretch_start ",stretch_start, " stretch_end ",stretch_end
	substring=sequence[stretch_start-1:stretch_end-1]
	o.write(">PUTATIVE_AA_"+str(stretch_start)+"_"+str(stretch_end)+"\n")
	o.write(substring+"\n")
#	maxi=max(zones)
#	centro=zones.index(maxi)

o.close()


