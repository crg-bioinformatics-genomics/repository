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
import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn import preprocessing
from sklearn import cross_validation
from sklearn import neighbors
from sklearn import metrics
from sklearn import svm
from sklearn import ensemble
from sklearn.externals import joblib



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


#def create_dic_fromFastaSingle(inputfilename):
#    infile=open(inputfilename,"r")
#    name=""
#    string=""
#    dic={}
#    for line in infile:	
#    	dic[line.split()[0]]=line.split()[1]
#    return dic

def create_dic_fromFasta(inputfilename):
    infile=open(inputfilename,"r")
    name=""
    string=""
    dic={}
    for line in infile:
	if line[0]==">":
		if (name==""):
			name=line.split()[0][1:]
		else:
			dic[name]=string			
			name=line.split()[0][1:]		
			string=""
	else:
		string+=line[:-1]
    dic[name]=string
    infile.close()
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
	
def bigFilter(RBP,scaleMeanStd):
	BIGmodel=joblib.load("data/MLmodels/FinalKNNModel.model")
	BIGmodel.classes_
	C025=C075=C050=0
	#print RBP
	data = numpy.zeros([1, 777])
	total=numpy.zeros([1, 259])
	#print data,total
	p=tp=0
	for key in (selectedMotSca):
		RRM=str(key.split()[0])
		scale=int(key.split()[1])
		errorthres=float(selectedMotSca[key].split()[3]) ####domain scales aucError bestThresError aucCorr bestCorr thresholdForAccum
		bestcorr=float(selectedMotSca[key].split()[5]) ####domain scales aucError bestThresError aucCorr bestCorr thresholdForAccum
		thresaccum=float(selectedMotSca[key].split()[6]) ####domain scales aucError bestThresError aucCorr bestCorr thresholdForAccum
		currentRRMprof=motifprofiles[RRM]
		currentRBPprof=posRBP_profiles[RBP][scale]
		protlen=len(posRBP[RBP])
		#calculate accumulation
		accumprofile= accum(currentRRMprof, currentRBPprof,protlen)
		if accumprofile!="":
			total[0,tp]=accumprofile.split()[int(bestcorr*10+10)]
			tp+=1
			C025=accumprofile.split()[int(0.25*10+10)] ##value at Correlation 0.25
			C075=accumprofile.split()[int(0.75*10+10)] ##value at Correlation 0.75
			C050=accumprofile.split()[int(0.50*10+10)] ##value at Correlation 0.50
			data[0,p]=float(C025)
			data[0,p+1]=float(C075)
			data[0,p+2]=float(C050)
			p+=3
		else:
			data[0,p]=0.0
			data[0,p+1]=0.0
			data[0,p+2]=0.0
			p+=3
			total[0,tp]=0.0
			tp+=1
	#print data,total
	prediction=BIGmodel.predict_proba(data)
	return total,prediction,BIGmodel.classes_
	
def smallFilter():
	smallModel=joblib.load('data/MLmodels/FinalRFModel.model')
	datafilter2 = numpy.zeros([1, 777])
	dicti={}
	position=10
	if (smallModel.classes_)[1]==1:
		position=1
	else:
		position=0
	infile=open("outputs/overlaps.txt","r") #RBPid RRMid PFAMid scale mystart myend maxcorr bestcorr minerr besterr
	count=0
	for line in infile:
		if line[0]!="#":
			Rid=line.split()[1]
			pf=line.split()[2]
			s=line.split()[3]
			maxcorr=line.split()[6]
			minerror=line.split()[8]
			mlVector= codify(Rid,s,maxcorr,error,selectedMotSca)
			#print mlVector
			datafilter2[0,]=numpy.array(mlVector)
			prediction=smallModel.predict_proba(datafilter2)
			ins=[]
			if prediction[0][position]>=0.5:
				ins.append(prediction[0][position])
				ins.append(line)
				dicti[count]=ins
				count+=1
	return dicti
			
	


def codify(RRM,scale,corr,err,MotSca): ##creates SVM vector
	piece1=[]
	piece2=[]
	piece3=[]
	for key in selectedMotSca:
	#	print key
		if RRM+" "+scale==key:
			piece1.append(1)
			piece2.append(corr)
			piece3.append(err)
		else:
			piece1.append(0)
			piece2.append(0)
			piece3.append(0)
	#print piece1
	#print piece2
	#print piece3
	return piece1+piece2+piece3


#######################main########################################

AAs = "ARNDCQEGHILKMFPSTWYVX"    
titles, property_table  = prof.parse_scales(  os.path.join("", "./scales/scales_generated.txt"), AAs )    


#####gest best motif+scales+correlation+error
selectedMotSca={}
ff=open("data/FinalError_Cand_LysateBilanciatoperCandNormalizedNEW/joinedCorrError.NEW","r") ###domain scales aucError bestThresError aucCorr bestCorr

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
#posRBP=create_dic_fromFastaSingle("data/RNAInKnown.fasta")
#posRBP=create_dic_fromFasta("data/lysate_30Cdhit.fasta")
#posRBP=create_dic_fromFastaSingle("data/inseq.fasta")
#posRBP=create_dic_fromFasta("data/candidate.cdhit30.fasta")
posRBP=create_dic_fromFasta("data/inseq.fasta")

posRBP_profiles={}

for key in posRBP:
    #print key,posRBP[key]
    smoothed, norm= cleverMachine_petr(posRBP[key],property_table,AAs)
    posRBP_profiles[key]=smoothed
print "Known RBPs done....",len(posRBP)," ",len(posRBP_profiles)

############calculate meand and stddev for each scale
scaleMeanStd={}
name=""
for i in range(len(titles)):
	scaleMeanStd[i]=str(numpy.mean(property_table[i][:-1]))+" "+str(numpy.std(property_table[i][:-1]))
logfile=open("log.log","w")

for RBP in posRBP:
################CHECK HMMER ##############################
	logfile.write(RBP+"\n")
	logfile.write("load pfamModels-pfamID mapping\n")
	infile=open("data/PfamOnlyRNAbinding.txt","r")
	pfamRNA={}
	for line in infile:#PDB_ID	CHAIN_ID	PdbResNumStart	PdbResNumEnd	PFAM_ACC	PFAM_Name	PFAM_desc	eValue
		pfamRNA[line.split()[4]]=line.split()[5] #model=pfamID
	#print pfam	

	o=open("data/hmmer.fasta","w")
	o.write(">"+str(RBP)+"\n")
	o.write(posRBP[RBP]+"\n")
	o.close()
	subprocess.call("hmmscan --cut_ga --noali --domtblout outputs/fragmentation.log /mnt/large/databases/Pfam_hmm/Pfam-A.hmm data/hmmer.fasta > outputs/hmm_log.txt",shell=True)
	infile=open("outputs/fragmentation.log","r")
	
	RRMfound=0
	for line in infile:
		if line[0]!="#":
		#	print line
			if line.split()[1] in pfamRNA:
				RRMfound+=1
			
	infile.close()
	##if HMMER does not find  RRMs (RRMfound<1) then proceed with the profiler 
	logfile.write("RRMfound "+str(RRMfound)+"\n")
	if RRMfound<1:	
	#	if True:
		totaldata = numpy.zeros([1, 259])
########################################FIRST filter BIGmodel########################################
		totaldata,prediction,classes= bigFilter(RBP,scaleMeanStd)
		print totaldata, prediction
		if (prediction[0][1]>prediction[0][0]):
			logfile.write("Passed\n")

		PfamID={}
		infile=open("data/pfam.id","r")#ENSP00000393530_136-212 IMLKRIYRSTPPEVIVEVLEPYVRLTTANVRIIKNRTGPMGHTYGFIDLDSHAEALRVVKILQNLDPPFSIDGKMVA PF14259.1
		for line in infile:
			PfamID[line.split()[0]]=line.split()[2]
		infile.close()
		scaleID={}
		infile=open("data/titlesofScales.txt","r")#scaleNumber scaleName
		for line in infile:
			scaleID[line.split()[0]]=line[:-1]
		infile.close()
		N=0
		
		errorCorr=[]
		keyposition=-1
		for key in (selectedMotSca):
			keyposition+=1
			RRM=str(key.split()[0])
			scale=int(key.split()[1])
			errorthres=float(selectedMotSca[key].split()[3]) ####domain scales aucError bestThresError aucCorr bestCorr thresholdForAccum
			bestcorr=float(selectedMotSca[key].split()[5]) ####domain scales aucError bestThresError aucCorr bestCorr thresholdForAccum
			thresaccum=float(selectedMotSca[key].split()[6]) ####domain scales aucError bestThresError aucCorr bestCorr thresholdForAccum
			if (bestcorr>=0.6): #!!YES o NO????			
				currentRRMprof=motifprofiles[RRM]#[scale]
				currentRBPprof=posRBP_profiles[RBP][scale]
				protlen=len(posRBP[RBP])
				RRMlen=len(motifs[RRM])
				RRMProfile=open("data/domainProfiles/"+RRM+"_"+str(scale)+".txt","w")
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
				#print totaldata,keyposition
				#accumprofile= accum(currentRRMprof, currentRBPprof,protlen)
				#if accumprofile!="":
				#	decisionvalue=accumprofile.split()[int(bestcorr*10+10)]
				if True:
		 			#if float(decisionvalue)< float(thresaccum):
					if totaldata[0,keyposition]<thresaccum:
						logfile.write(RRM+" "+str(scale)+" negative\n")
						N+=1
					else:
						logfile.write(RRM+" "+str(scale)+" positive\n")
						#if (len(currentRRMprof) < len(currentRBPprof)):
						shifts=[]
						for sw in range(0,len(currentRBPprof)-len(currentRRMprof)):
							####calculate Error####################
							error=0
							y=normalize(currentRRMprof,scale,scaleMeanStd)
							x=normalize(currentRBPprof[sw:sw+len(currentRRMprof)],scale,scaleMeanStd)
							error= numpy.sqrt(numpy.sum((x-y)**2))/RRMlen
							####calculate Correlation####################
							correlation=float(numpy.corrcoef(currentRRMprof,currentRBPprof[sw:sw+len(currentRRMprof)])[0, 1])
							if (correlation>=float(bestcorr)) and (float(error)<=errorthres) and (correlation>=0.7):
								#runvar is key
								shifts.append(RBP+" "+RRM+" "+str(scale)+" "+str(sw)+" "+str(correlation)+" "+str(bestcorr)+" "+str(error)+" "+str(errorthres))		
						swold=0
						trip=[]
						for i in shifts:
							if trip==[]:

								trip.append(i)
								swold=int(i.split()[3])
							else:
								swnow=int(i.split()[3])
								if (swold+1)==swnow:
									trip.append(i)
									swold=swnow
								else:
									if len(trip)>=3:
										for xy in trip:
											errorCorr.append(xy)
										trip=[]
										trip.append(i)
										swold=int(i.split()[3])
						
									else:
										trip=[]
										trip.append(i)
										swold=int(i.split()[3])	
						if len(trip)>=3:
							for xy in trip:
								errorCorr.append(xy)
				else:
					logfile.write("RRM longer than protein\n")
			else:
				logfile.write(RRM+" correlation lower 0.6\n")
				#N+=1
		logfile.write("create error_corr file\n")
		error_corr=open("outputs/errorCorrRBP_RRM.txt","w")
		error_corr.write("#RBP RRM PfamID scale sw corr corrthres error errorthres\n")		
		#print errorCorr
		countimmage=1
		#rbpout={}
		for a in errorCorr:
			#entry=[]
			for i in range(len(a.split())):	
				if i==1:
					error_corr.write(a.split()[i]+" "+PfamID[a.split()[i]]+" ")
					#if (PfamID[a.split()[1]]+" "+scaleID[a.split()[2]]) not in rbpout:
					#	countimmage+=1
					#	entry.append(scaleID[a.split()[2]])
					#	entry.append(PfamID[a.split()[1]])
					#	entry.append(a.split()[3])
					#	entry.append(a.split()[4])
					#	rbpout[a.split()[1]+" "+a.split()[2]]=entry#+" Image"+str(countimmage)
					#else:
					#	entry.append(scaleID[a.split()[2]])
					#	entry.append(PfamID[a.split()[1]])
					#	entry.append(a.split()[3])
					#	entry.append(a.split()[4])
					#	test=rbpout[a.split()[1]+" "+a.split()[2]]
					#	if float(test.split()[1])< float(a.split()[4]):
			#				rbpout[a.split()[1]+" "+a.split()[2]]=entry
				else:	
					error_corr.write(a.split()[i]+" ")
			error_corr.write("\n")		
		error_corr.close()	
		if (sum(1 for line in open('outputs/errorCorrRBP_RRM.txt')))>1:##if at least one profile has been detected!!
			rbpout={}	
			logfile.write("second filter \n")
			#######SECOND Filter################	
			###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CONTINUE HERE!!!!!!!!!!!!!!!!!!!
			subprocess.call("cat CheckOverlap2.R | R --slave --vanilla ",shell=True)
			prediction=smallFilter()
			if prediction!={}:
				o=open("outputs/errorCorrRBP_RRM.Final","w")
				for key in prediction:
						o.write(prediction[key][1])	
				o.close()
				logfile.write("second filter... OK\n")
				###################create plot output########################################################
				logfile.write("create plots\n")
				#subprocess.call("cat GetResultsstraightplotsBeta.R | R --slave --vanilla --args "+str(protlen),shell=True)		
				#subprocess.call("cat Rscript.R | R --slave --vanilla",shell=True)			
				subprocess.call("cat Plots.R | R --slave --vanilla --args "+str(protlen),shell=True)		
				logfile.write("create plots OK\n")
				###################create HTML table for output########################################################
				tablefile=open("outputs/testtable.html","w")
				tablefile.write("<tbody>\n")
				for key in prediction: ##creates html table with: region-length pfamID scale corr
					#print key, rbpout[key]
					tablefile.write("<tr>\n")	
					RRMlength=sum(1 for line in open("data/domainProfiles/"+prediction[key][1].split()[1]+"_"+str(prediction[key][1].split()[3])+".txt"))
					tablefile.write("<td>"+str(prediction[key][1].split()[4])+"-"+str(int(prediction[key][1].split()[5])+RRMlength)+"</td>\n") ##region					
					tablefile.write("<td>"+str(prediction[key][1].split()[2])+"</td>\n") ##pfam
					tablefile.write("<td>"+str(round(float(prediction[key][1].split()[6]),2))+"</td>\n") ##correlation
					tablefile.write("<td>"+scaleID[prediction[key][1].split()[3]]+"</td>\n") #scalename
					tablefile.write("<td> <image src=\""+str(prediction[key][1].split()[1])+"_"+str(prediction[key][1].split()[3])+".jpeg\" width=\"80%\" height=\"20%\" /> </td>\n")
					tablefile.write("<td> <h3><a href=\""+str(prediction[key][1].split()[1])+"_"+str(prediction[key][1].split()[3])+".jpeg\">Download jpeg</a></h3>  </td></tr>\n")
				tablefile.write("</tbody>\n")
				tablefile.close()
			else:
				#subprocess.call("awk \'($1!=\"#\"&&$3!=\"-----\"){print $2\" \"$3\" \"$4\" \"$6\" \"$18\" \"$19}\' outputs/fragmentation.log >outputs/preparedFrag.log",shell=True)
				logfile.write(RBP+" No domain found\n")	
				subprocess.call("cat NoDomainFound.R | R --slave --vanilla ",shell=True)
				logfile.write(RBP+" No domain found OK\n")	


			#subprocess.call("awk \'($1!=\"#\"&&$3!=\"-----\"){print $2\" \"$3\" \"$4\" \"$6\" \"$18\" \"$19}\' outputs/fragmentation.log >outputs/preparedFrag.log",shell=True)
			#subprocess.call("cat NoPass.R | R --slave --vanilla ",shell=True)
			#subprocess.call("cat CheckOverlap.R | R --slave --vanilla ",shell=True)
		else:
			#subprocess.call("awk \'($1!=\"#\"&&$3!=\"-----\"){print $2\" \"$3\" \"$4\" \"$6\" \"$18\" \"$19}\' outputs/fragmentation.log >outputs/preparedFrag.log",shell=True)
			logfile.write(RBP+" No domain found\n")	
			subprocess.call("cat NoDomainFound.R | R --slave --vanilla ",shell=True)
			logfile.write(RBP+" No domain found OK\n")	
	else:
		##if HMMER finds  RRMs stop the profiler and plot the detected domains
		logfile.write(RBP+" no pass\n")
		subprocess.call("awk \'($1!=\"#\"&&$3!=\"-----\"){print $2\" \"$3\" \"$4\" \"$6\" \"$18\" \"$19}\' outputs/fragmentation.log >outputs/preparedFrag.log",shell=True)
		logfile.write(RBP+" no pass...1 OK\n")
		subprocess.call("cat NoPass.R | R --slave --vanilla ",shell=True)
		logfile.write(RBP+" no pass.. 2 OK\n")
