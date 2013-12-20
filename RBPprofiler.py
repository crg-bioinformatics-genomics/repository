from __future__ import division
import os
import sys
import random
import subprocess
import pickle
import math
import numpy
import scipy
import collections
from matplotlib import pyplot as plt
# include CM modules
import libprofiles as prof
import libcombinations as lc  


def cleverMachine_petr(test_sequence,property_table,AAs):
    """
    profile calculations from Petr 
    """
    def process_sequence(sequence, window_size=7):
        """Single sequence evaluation using scales provided"""
        records = [("seq_submit", sequence)]  # + sequences
        scored = prof.score_proteins(records, AAs, property_table)
        smoothed = map(lambda scored: prof.smooth(scored, window_size=window_size), scored)
        cur_prot = smoothed[0]
        return cur_prot,scored
    # derives scale information
    NO_OF_PROPERTIES = property_table.shape[0]
    # property_groups, ungrouped_properties = lc.property_group_init(NO_OF_PROPERTIES)
    # contains 80x379 entries - 80 scales with 379 entries in each row
    result_smoothed ,prob_resultsnorm= process_sequence(test_sequence)
    return result_smoothed,prob_resultsnorm


def create_dic_fromFasta(inputfilename):
    #log.write(inputfilename)
    #print ("Reading sequences...")
    infile=open(inputfilename,"r")
    name=""
    string=""
    start=0
    dic={}
    num=0
    for line in infile:
        #print line
        if line[0]==">":
            if start==1:
                if string!="Sequence unavailable" and string!="No UUR is annotated for this transcript":
                    dic[name]=string
                    num+=1
                #print name, num
                string=""
                name=""
                new=""
                name=line[0:len(line)]

                for a in name:
                    if a!=">" and a!="\n":
                        new+=a
                        #let=raw_input(": ")
                name=new
            else:
                start=1
                new=""
                name=line[0:len(line)]
                for a in name:
                    if a!=">" and a!="\n":
                        new+=a
                        #let=raw_input(": ")
                name=new
        else:
            for a in line:
                if a!="\n":
                    if a=="T":
                        string+="U"
                    else:
                        string+=a

    if string!="Sequence unavailable" and string!="No UUR is annotated for this transcript":
        dic[name]=string
        num+=1
  
    return dic
    


AAs = "ARNDCQEGHILKMFPSTWYVX"    
titles, property_table  = prof.parse_scales(  os.path.join("", "./scales/scales_generated.txt"), AAs )    
#print property_table

log=open("./output/log.log","w")

#####gest best motif+scales
selectedMotSca={}
#ff=open("tmpu3.txt","r")
ff=open("./data/selectedScales.txt","r")
for line in ff:
	#print line
	selectedMotSca[line]=(line.split("+")[1][:-1])+" "+(line.split("+")[0])
#	print(line.split("+")[1][:-1])+" "+(line.split("+")[0])
log.write(str(len(selectedMotSca))+"\n")
#let=raw_input(": ")

######calculate motif profile
motifs={}
motifs=create_dic_fromFasta("./data/motifs.txt")
motifprofiles={}
for key in motifs:
	for a in selectedMotSca:
		mm=a.split("+")[0]
		ss=int(a.split("+")[1])
		if mm==key:
	    #print "yes",key
		    smoothed, norm= cleverMachine_petr(motifs[key],property_table,AAs)
		    motifprofiles[key]=smoothed[ss]
log.write(str(len(motifprofiles))+" from Motifs "+str(len(motifs))+" selected\n")
#let=raw_input(": ")    
    
######calculate candidate profile
#unknownRBP={}
#unknownRBP=create_dic_fromFasta("candidate.cdhit30.fasta")
#unknownRBP_profiles={}
#for key in unknownRBP:
#    smoothed, norm= cleverMachine_petr(unknownRBP[key],property_table,AAs)
#    unknownRBP_profiles[key]=smoothed
#print "Candidates done....", len(unknownRBP),len(unknownRBP_profiles)

######calculate known RBP profile
posRBP={}
posRBP=create_dic_fromFasta("./data/inseq.fasta")
posRBP_profiles={}

for key in posRBP:
    #print key,posRBP[key]
    smoothed, norm= cleverMachine_petr(posRBP[key],property_table,AAs)
    posRBP_profiles[key]=smoothed
#print posRBP_profiles
log.write("Known RBPs done...."+str(len(posRBP))+" "+str(len(posRBP_profiles))+"\n")

outfile=open("./output/profileForProtein.txt","w")
for RBP in posRBP:
	log.write(RBP+"\n")
	for key in (selectedMotSca):
#		print key
		
		RRM=str(key.split("+")[0])
		scale=int(key.split("+")[1])
		#print RRM,scale
		#let=raw_input(": ")
		
		log.write(RRM+"\t"+str(scale)+"\n")
		#print RRM,scale
		
		currentRRMprof=motifprofiles[RRM]#[scale]
		currentRBPprof=posRBP_profiles[RBP][scale]
		protlen=len(posRBP[RBP])
		RRMlen=len(motifs[RRM])
		#          
		finalcorr=[] 
		correlationField=[]
		if (len(currentRRMprof))<(len(currentRBPprof)):
			for sw in range(0,len(currentRBPprof)-len(currentRRMprof)):
				correlation=numpy.corrcoef(currentRRMprof,currentRBPprof[sw:sw+len(currentRRMprof)] )[0, 1]
				try:
					correlationField.append(round(correlation,2))
				except:
					correlationField=round(numpy.nan_to_num(correlation),2)    
			log.write(RRM+str(scale)+"\n")
			#print correlationField
			for i in range(-10,11):
				icorr=i/10
				#print icorr
				counts=0
				for a in correlationField: 
					if a<icorr:
						counts+=1
				finalcorr.append(counts)	

			#print finalcorr
			#let=raw_input(": ")
			outstring=RBP+"\t"+RRM+"\t"+str(scale)+"\t"
			for ii in finalcorr:
				outstring+=str(round(int(ii)/int(protlen),2))+"\t"
			outstring+=str(protlen)+"\t"+str(RRMlen)+"\n"
			#print outstring
			outfile.write(outstring)
			#print outstring
outfile.close()

#let=raw_input(": ")
log.write("execute Rscript\n")
subprocess.call("cat Rscript.R | R --slave --vanilla",shell=True)
log.write("done\n")

log.write("create profile...\n")
pfam=open("./data/pfam.id","r")#ENSP00000393530_136-212 IMLKRIYRSTPPEVIVEVLEPYVRLTTANVRIIKNRTGPMGHTYGFIDLDSHAEALRVVKILQNLDPPFSIDGKMVA PF14259.1
resultfile=open("./output/RBProfile.final","r") # "ENSP00000437886_261-336 17 0.936"
newmotif="error"
motif="error"
scale="error"
score="error"
for line in resultfile:
	motif=line.split()[0]
	scale=int(line.split()[1])
	score=float(line.split()[2])
	for pf in pfam:
		if pf.split()[0]==motif:
			newmotif=pf.split()[2]
			log.write("Pfam substituted\n")
pfam.close()
resultfile.close()
#print motif,scale,score,newmotif

currentRRMprof=motifprofiles[RRM]
for RBP in posRBP:
	currentRBPprof=posRBP_profiles[RBP][scale]
	finalcorr=[] 
	correlationField=[]
	if (len(currentRRMprof))<(len(currentRBPprof)):
		for sw in range(0,len(currentRBPprof)-len(currentRRMprof)):
			correlation=numpy.corrcoef(currentRRMprof,currentRBPprof[sw:sw+len(currentRRMprof)] )[0, 1]
			try:
				correlationField.append(round(correlation,2))
			except:
				correlationField=round(numpy.nan_to_num(correlation),2)    
outstring=""
for a in correlationField:
	outstring+=str(a)+" "
for a in range(len(outstring.split()),protlen):
	outstring+=str(0)+" "

outfile=open("./output/RBP.out","w")
outfile.write(outstring+"\n")
#print protlen,len(outstring.split()),RRMlen
#let=raw_input(": ")
outfile.write("#domain: "+newmotif+"\n")
outfile.write("#scale: "+titles[scale]+"\n")
outfile.write("#score: "+str(score)+"\n")
outfile.close()

log.write("create Immage...")
subprocess.call("cat createImmage.R | R --slave --vanilla --args "+str(score)+" "+newmotif,shell=True)
		
log.close()


