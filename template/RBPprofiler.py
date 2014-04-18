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
                	string+=a

    if string!="Sequence unavailable" and string!="No UUR is annotated for this transcript":
        dic[name]=string
        num+=1
  
    return dic


def create_dic_fromFastaSingle(inputfilename):
    infile=open(inputfilename,"r")
    name=""
    string=""
    dic={}
    for line in infile:
        if line[0]==">":
        	if len(line.split())>1:
			name=line.split()[0][:-1]
			string=""
		else:
			name=line[1:len(line)-1]
			string=""
        else:
	    if "\n" in line:
	        string+=line[:-1]
	    else:
		string+=line
            	
    	dic[name]=string
    return dic


def create_dic_fromFastaOneline(inputfilename):
    infile=open(inputfilename,"r")
    name=""
    string=""
    dic={}
    for line in infile:
	name=line.split()[0]
	stringi=line.split()[1]
	dic[name]=stringi
    return dic


def normalize(profile,scale,scaleMeanStd):
##zeta normalization
	result=scaleMeanStd[scale]
	mean=float(result.split()[0])
	std=float(result.split()[1])
	new=[]
	for a in profile:
		new.append(round((a-mean)/std,2))
	return 	numpy.array(new)


def accum(currentRRMprof, currentRBPprof,protLen):
##calculates the accum of the profile
	correlationField={}
	for a in range(21): #-1...0...1):
		correlationField[a]=0
	outstring=""
	#print correlationField, currentRRMprof, currentRBPprof
	for key in correlationField: ##init dict
		correlationField[key]=0
	for sw in range(0,len(currentRBPprof)-len(currentRRMprof)):
		correlation=numpy.corrcoef(currentRRMprof,currentRBPprof[sw:sw+len(currentRRMprof)] )[0, 1]
		try:
			enterinField=round(correlation,1)
			correlationField[abs(enterinField*10+10)]+=1
		except:
			enterinField=round(numpy.nan_to_num(correlation),1)
			correlationField[(enterinField*10+10)]+=1  	    
		od = collections.OrderedDict(sorted(correlationField.items()))
		dictentry=""
 		for k, v in od.iteritems(): 
			dictentry+=str(v)+"\t"
		accumula=int(dictentry.split()[0])
		outstring=""
		#print accumula, len(dictentry.split()),protLen
		for i in range(len(dictentry.split())):
			accumula+=round((int(dictentry.split()[i])/int(protLen)),2)
	#		accumula+=float(dictentry.split()[i])
			outstring+=str(accumula)+" "
	return outstring
		
	


#######################main########################################

AAs = "ARNDCQEGHILKMFPSTWYVX"    
titles, property_table  = prof.parse_scales(  os.path.join("", "./scales/scales_generated.txt"), AAs )    


#####gest best motif+scales+correlation+error
selectedMotSca={}
ff=open("data/FinalError_Cand_LysateBilanciatoperCand/joinedCorrError.txt","r") ###domain scales aucError bestThresError aucCorr bestCorr

for line in ff:
	if line[0]!="#":
		motif=line.split()[0]
		scale=line.split()[1]
		selectedMotSca[motif+" "+scale]=line[:-1]
print len(selectedMotSca)


######calculate motif profile for preselected domain+scale
motifs={}
motifs=create_dic_fromFasta("data/motifs.txt")
motifprofiles={}
for key in motifs:
	for a in selectedMotSca:
		mm=a.split()[0]
		ss=int(a.split()[1])
		if mm==key:
		    smoothed, norm= cleverMachine_petr(motifs[key],property_table,AAs)
		    motifprofiles[key]=smoothed[ss]

print len(motifprofiles)," from Motifs ",len(motifs)," selected"
#print motifs

######calculate known RBP profile
posRBP={}
#posRBP=create_dic_fromFasta("../rawData/Gutfreund_negative.fasta")
#posRBP=create_dic_fromFasta("data/known.cdhit30.fasta")
posRBP=create_dic_fromFastaOneline("data/RNAInKnown.fasta")
print posRBP
#posRBP=create_dic_fromFastaSingle("data/lysate_30Cdhit.fasta")
#posRBP=create_dic_fromFastaSingle("data/inseq.fasta")
#posRBP=create_dic_fromFasta("data/candidate.cdhit30.fasta")
#posRBP=create_dic_fromFasta("inseqLysate.neg")
#posRBP=create_dic_fromFasta("inseqCandidate.fasta")
posRBP_profiles={}

for key in posRBP:
    #print key,posRBP[key]
    smoothed, norm= cleverMachine_petr(posRBP[key],property_table,AAs)
    posRBP_profiles[key]=smoothed
print "Known RBPs done....",len(posRBP)," ",len(posRBP_profiles)
############calculate meand and stddev for each scale

#ENSP00000366843

scaleMeanStd={}
name=""
for i in range(len(titles)):
	scaleMeanStd[i]=str(numpy.mean(property_table[i][:-1]))+" "+str(numpy.std(property_table[i][:-1]))

for RBP in posRBP:
	print RBP


	print "load pfamModels-pfamID mapping"
	infile=open("data/PfamOnlyRNAbinding.txt","r")
	pfam={}
	for line in infile:#PDB_ID	CHAIN_ID	PdbResNumStart	PdbResNumEnd	PFAM_ACC	PFAM_Name	PFAM_desc	eValue
		pfam[line.split()[4]]=line.split()[5] #model=pfamID
	#print pfam	

	o=open("data/input_sequence.fasta","w")
	o.write(">"+str(RBP)+"\n")
	o.write(posRBP[RBP]+"\n")
	o.close()
	subprocess.call("hmmscan --cut_ga --noali --domtblout outputs/fragmentation.log /mnt/large/databases/Pfam_hmm/Pfam-A.hmm data/input_sequence.fasta > outputs/hmm_log.txt",shell=True)
	infile=open("outputs/fragmentation.log","r")
	
	RRMfound=0
	for line in infile:
		if line[0]!="#":
		#	print line
			if line.split()[1] in pfam:
				RRMfound+=1
			
	infile.close()
	
	print "RRMfound",RRMfound
#	if RRMfound<1:	
	if True:
		#outfile=open("data/SVMtest.te","w")
		##CHECK filter SVM
		#index=0
		#outstring=""
		#for key in (selectedMotSca):
		##	print key
		#	index+=1	
		#	RRM=str(key.split()[0])
		#	scale=int(key.split()[1])
		#	errorthres=float(selectedMotSca[key].split()[3]) ####domain scales aucError bestThresError aucCorr bestCorr thresholdForAccum
		#	bestcorr=float(selectedMotSca[key].split()[5]) ####domain scales aucError bestThresError aucCorr bestCorr thresholdForAccum
		#	thresaccum=float(selectedMotSca[key].split()[6]) ####domain scales aucError bestThresError aucCorr bestCorr thresholdForAccum
		#	currentRRMprof=motifprofiles[RRM]#[scale]
		#	currentRBPprof=posRBP_profiles[RBP][scale]
		#	protlen=len(posRBP[RBP])
		#	####calculate accum ###################	
		#	accumprofile= accum(currentRRMprof, currentRBPprof,protlen)
		#	if accumprofile!="":
		#		#print accumprofile
		#		decisionvalue=accumprofile.split()[int(bestcorr*10+10)]
		#		#print "value at correlation",bestcorr," is ",decisionvalue,". decision is at", thresaccum
		#		outstring+=str(index)+":"+str(decisionvalue)+" "
		#if len(outstring)>2:
		#	outfile.write(outstring+"\n")
		#outfile.close()
		#subprocess.call("./data/svm-scale -r data/fileFinal  data/SVMtest.te > data/SVMscaled.te",shell=True)
		#subprocess.call("./data/svm-predict -b 1 -q data/SVMscaled.te data/Finalmodel.m data/test.out",shell=True)
		#ii=open("data/test.out","r")
		#for line in ii:
		#	if line.split()[0]!="labels":
		#		go=int(line.split()[0])

		PfamID={}
		infile=open("data/pfam.id","r")#ENSP00000393530_136-212 IMLKRIYRSTPPEVIVEVLEPYVRLTTANVRIIKNRTGPMGHTYGFIDLDSHAEALRVVKILQNLDPPFSIDGKMVA PF14259.1
		for line in infile:
			PfamID[line.split()[0]]=line.split()[2]
		infile.close()
		scaleID={}
		infile=open("data/titlesofScales.txt","r")#ENSP00000393530_136-212 IMLKRIYRSTPPEVIVEVLEPYVRLTTANVRIIKNRTGPMGHTYGFIDLDSHAEALRVVKILQNLDPPFSIDGKMVA PF14259.1
		for line in infile:
			scaleID[line.split()[0]]=line[:-1]
		infile.close()
		print "passed"
		N=0
		
		errorCorr=[]
		for key in (selectedMotSca):
	#		print key
		
			RRM=str(key.split()[0])
			scale=int(key.split()[1])
			errorthres=float(selectedMotSca[key].split()[3]) ####domain scales aucError bestThresError aucCorr bestCorr thresholdForAccum
			bestcorr=float(selectedMotSca[key].split()[5]) ####domain scales aucError bestThresError aucCorr bestCorr thresholdForAccum
			thresaccum=float(selectedMotSca[key].split()[6]) ####domain scales aucError bestThresError aucCorr bestCorr thresholdForAccum
			#print RBP,RRM,scale, errorthres,bestcorr, thresaccum
			#let=raw_input(": ")
			if bestcorr>=0.6:
				#print RRM,str(scale)
				#print RRM,scale
		
				currentRRMprof=motifprofiles[RRM]#[scale]
				currentRBPprof=posRBP_profiles[RBP][scale]
				protlen=len(posRBP[RBP])
				RRMlen=len(motifs[RRM])
				RRMProfile=open("data/domainProfiles/"+str(PfamID[RRM])+"_"+str(scale)+".txt","w")
				#for i in normalize(currentRRMprof,scale,scaleMeanStd):
				h=0
				for i in currentRRMprof:
					if h==0	:
						for t in range(3):			
							RRMProfile.write(str(i)+"\n")			
						h=1
					RRMProfile.write(str(i)+"\n")			
				for t in range(3):
					RRMProfile.write(str(i)+"\n")			
				RRMProfile.close()

				#	subprocess.call("R CMD BATCH Rscript.R",shell=True)
				#let=raw_input(": ")

				#print RRM, scale, errorthres,corrthres
				proteinProfile=open("data/proteinProfiles/"+RBP+"_"+str(scale)+".txt","w")
				h=0
				for i in currentRBPprof:
				#for i in normalize(currentRBPprof,scale,scaleMeanStd):
					if h==0:
						for t in range(3):
							proteinProfile.write(str(i)+"\n")
						h=1		
					proteinProfile.write(str(i)+"\n")		
				for t in range(3):
					proteinProfile.write(str(i)+"\n")
				proteinProfile.close()
				###calculate accum ###################
		
				accumprofile= accum(currentRRMprof, currentRBPprof,protlen)
				if accumprofile!="":
					#print accumprofile
					decisionvalue=accumprofile.split()[int(bestcorr*10+10)]
					#print "value at correlation",bestcorr," is ",decisionvalue,". decision is at", thresaccum
		 			if float(decisionvalue)< float(thresaccum):
						print RRM,scale," negative"
						N+=1
					else:
						print RRM,scale,"  positive"
						#if (len(currentRRMprof) < len(currentRBPprof)):
						shifts=[]
						for sw in range(0,len(currentRBPprof)-len(currentRRMprof)):
						####calculate Error####################
							error=0
							y=normalize(currentRRMprof,scale,scaleMeanStd)
							x=normalize(currentRBPprof[sw:sw+len(currentRRMprof)],scale,scaleMeanStd)
							error= numpy.sqrt(numpy.sum((x-y)**2))/RRMlen
						#	let=raw_input(": ")
			
							####calculate Correlation####################
							correlation=float(numpy.corrcoef(currentRRMprof,currentRBPprof[sw:sw+len(currentRRMprof)])[0, 1])
							if (correlation>=float(bestcorr)) and (float(error)<=errorthres) and (correlation>=0.7):
								#runvar is key
								shifts.append(RBP+" "+RRM+" "+str(scale)+" "+str(sw)+" "+str(correlation)+" "+str(bestcorr)+" "+str(error)+" "+str(errorthres))
	#					for i in shifts:
	#						print "orig",i
	#					let=raw_input(": ")					
						swold=0
						trip=[]
						#print shifts
						for i in shifts:
							#print "i",i
							if trip==[]:

								trip.append(i)
								swold=int(i.split()[3])
								#print "new", swold
							else:
								swnow=int(i.split()[3])
								#print "extend with",swnow
								if (swold+1)==swnow:

								#	print"ok"
									trip.append(i)
									swold=swnow
								else:
									#print trip
								#	print "finish"
								
									if len(trip)>=3:
										#print "write in file"
										for xy in trip:
											#error_corr.write(xy+"\n")
											errorCorr.append(xy)
										trip=[]
										trip.append(i)
										swold=int(i.split()[3])
						
									else:
										#print "nothing"
										#for xy in trip:
										#	print xy
										trip=[]
										trip.append(i)
										swold=int(i.split()[3])	
						
									#let=raw_input(": ")
					#	print "remaining",trip
						if len(trip)>=3:
							for xy in trip:
								#error_corr.write(xy+"\n")					
								errorCorr.append(xy)
					

				else:
					print "RRM longer than protein"
	#		else:
	#			print RRM+" correlation lower 0.6"
				#N+=1
		error_corr=open("outputs/errorCorrRBP_RRM.txt","w")
		error_corr.write("#RBP RRM scale sw corr corrthres error errorthres\n")		
		
		#print errorCorr
		countimmage=1
		rbpout={}

		for a in errorCorr:
			for i in range(len(a.split())):	
				if i==1:
					error_corr.write(PfamID[a.split()[i]]+" ")
					if (PfamID[a.split()[1]]+" "+scaleID[a.split()[2]]) not in rbpout:
						countimmage+=1
						rbpout[PfamID[a.split()[1]]+" "+scaleID[a.split()[2]]]=a.split()[3]+" "+a.split()[4]+" "+a.split()[2]#+" Image"+str(countimmage)
					else:
						entry=rbpout[PfamID[a.split()[1]]+" "+scaleID[a.split()[2]]]
						if float(entry.split()[1])< float(a.split()[4]):
							rbpout[PfamID[a.split()[1]]+" "+scaleID[a.split()[2]]]=a.split()[3]+" "+a.split()[4]+" "+a.split()[2]
							
				else:	
					error_corr.write(a.split()[i]+" ")
			error_corr.write("\n")

		
		error_corr.close()	
		subprocess.call("cat Rscript.R | R --slave --vanilla",shell=True)
		
	#	let=raw_input(": ")
###################create HTML table for output########################################################
		tablefile=open("outputs/testtable.html","w")
		tablefile.write("<tbody>\n")
		for key in rbpout: ##creates html table with: region-length pfamID scale corr
			#print key, rbpout[key]
			tablefile.write("<tr>\n")	
			tablefile.write("<td>"+str(rbpout[key].split()[0])+"</td>\n") 					
			tablefile.write("<td>"+str(key.split()[0])+"</td>\n") ##pfam
			rest=""
			for i in range(2,len(key.split())):			
				#print i
				rest+=key.split()[i]+" "  #scale
			tablefile.write("<td>"+str(rest)+"</td>\n")
			tablefile.write("<td>"+str(round(float(rbpout[key].split()[1]),2))+"</td>\n") 
			# <h3><a href="profiles.jpeg">Download jpeg</a></h3>
#			tablefile.write("<td> <image src=\""+str(rbpout[key].split()[2])+".jpeg\" /> </td>\n")# <image src="profiles.jpeg" /></br></br>
#			tablefile.write("<td> <h3><a href=\""+str(rbpout[key].split()[2])+".jpeg\">Download jpeg</a></h3>  </td></tr>\n")
			#<img src="img.jpg" width="50%" height="50%"/>
			tablefile.write("<td> <image src=\""+str(key.split()[0])+"_"+str(rbpout[key].split()[2])+".jpeg\" width=\"80%\" height=\"20%\" /> </td>\n")# <image src="profiles.jpeg" /></br></br>
#			tablefile.write("<td> <image src=\""+str(key.split()[0])+".jpeg\" /> </td>\n")# <image src="profiles.jpeg" /></br></br>
			tablefile.write("<td> <h3><a href=\""+str(key.split()[0])+"_"+str(rbpout[key].split()[2])+".jpeg\">Download jpeg</a></h3>  </td></tr>\n")

		tablefile.write("</tbody>\n")
		tablefile.close()
			
#		let=raw_input(": ")
		#print RBP, protlen,N,len(pfield)
###################create plot output########################################################
		subprocess.call("cat GetResultsstraightplots.R | R --slave --vanilla --args "+str(protlen),shell=True)	
			
		subprocess.call("awk \'($1!=\"#\"&&$3!=\"-----\"){print $2\" \"$3\" \"$4\" \"$6\" \"$18\" \"$19}\' outputs/fragmentation.log >outputs/preparedFrag.log",shell=True)
		subprocess.call("cat NoPass.R | R --slave --vanilla ",shell=True)
		subprocess.call("cat CheckOverlap.R | R --slave --vanilla ",shell=True)
		

#		let=raw_input(": ")

	else:
		subprocess.call("awk \'($1!=\"#\"&&$3!=\"-----\"){print $2\" \"$3\" \"$4\" \"$6\" \"$18\" \"$19}\' outputs/fragmentation.log >outputs/preparedFrag.log",shell=True)
		subprocess.call("cat NoPass.R | R --slave --vanilla ",shell=True)
		
		print RBP," no pass"
