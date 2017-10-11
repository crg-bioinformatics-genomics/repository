from __future__ import division
import os
import sys

sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/local/lib/python2.7/site-packages/pip-1.1-py2.7.egg')
sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/site-packages/pip-1.1-py2.7.egg')
sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7')
sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/plat-linux2')
sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/lib-tk')
sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/lib-old')
sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/lib-dynload')
sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/local/lib/python2.7/site-packages')
sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/site-packages')
sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/local/lib/python2.7/site-packages/IPython/extensions')



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
import smtplib


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
                        string+="T"
                    else:
                        string+=a

    if string!="Sequence unavailable" and string!="No UUR is annotated for this transcript":
        dic[name]=string
        num+=1
  
    return dic
 	

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
		liste=[]
		for i in range(21):
			accumula+=float(dictentry.split()[i])
			liste.append(accumula/int(protLen))
		newliste=[]
		for i in range(len(liste)):
			if (liste[19]!=0):
				newliste.append(liste[i]*100/(liste[19]*100))
			else:
				newliste.append(0)
		
		for i in newliste:
			outstring+=str(i)+" "

	return outstring
   
def features(RBP,RBPlen,RBP_profiles,scalet,motifprofiles):

	onlyBestCorr=numpy.zeros([1, len(scalet)])
	p=tp=0	
	for key in (scalet):
		RRM=str(key.split()[0])
		scale=int(key.split()[1])
		bestcorr=float(scalet[key].split()[1])
		#print key,RRM,scale,bestcorr,RBP
		rbpprof=RBP_profiles[RBP][scale]
		#print RBP_profiles[RBP][scale]
		#print motifprofiles
		rrmprof=motifprofiles[RRM+" "+str(scale)]
		
		#calculate accumulation
		accumprofile= accum(rrmprof, rbpprof,RBPlen)
		if accumprofile!="":
			onlyBestCorr[0,tp]=accumprofile.split()[int(bestcorr*10+10)]
			tp+=1
		else:
			onlyBestCorr[0,tp]=0.0
	#print data,onlyBestCorr
	#prediction=BIGmodel.predict_proba(data)
	return onlyBestCorr#,prediction,BIGmodel.classes_



AAs = "ARNDCQEGHILKMFPSTWYVX"    
titles, property_table  = prof.parse_scales(  os.path.join("", "./scales/scales_generated.txt"), AAs )    
#print titles

print "load scan files"
print "load C_scales"
selectedMotSca_short={}
infilename = "./data/scaletta_scan/C_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca_short["C "+motif+" "+scale]=corr

print "load NC_scales"
infilename = "./data/scaletta_scan/NC_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca_short["NC "+motif+" "+scale]=corr

print "load P_scales"
infilename ="./data/scaletta_scan/P_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca_short["P "+motif+" "+scale]=corr


print "load meta files"
print "load C_scales"
selectedMotSca_long={}
infilename = "./data/scaletta_meta/C_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca_long["C "+motif+" "+scale]=corr

print "load NC_scales"
infilename = "./data/scaletta_meta/NC_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca_long["NC "+motif+" "+scale]=corr

print "load P_scales"
infilename ="./data/scaletta_meta/P_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca_long["P "+motif+" "+scale]=corr



print "load pfamModels-pfamID mapping"
infile=open("./data/PfamIDsRnaBinding.txt","r")
cpfam={}
for line in infile:#PFAMID=1
	cpfam[line.split()[0].split(".")[0]]=1


print "load pfamModels-nonclassical mapping"
infile=open("./data/listOfNonClassicalBaltzCastelloKwon.Pfam","r")
ncpfam={}
for line in infile:#PFAMID = 1
	ncpfam[line.split()[0]]=1


print "load pfamModels-putative mapping"
infile=open("./data/listOfunknownPutativeBaltzCastelloKwon.Pfam","r")
putpfam={}
for line in infile:#PFAMID = 1
	putpfam[line.split()[0]]=1

print "load additional-pfamID mapping"
infile=open("./data/additional.Pfam","r")
addpfam={}
for line in infile:#PFAMID=1
	addpfam[line.split()[0]]=1



print "load motifs scan"
motifsshort={}
motifsshort=create_dic_fromFasta("./data/motifs.txt")
motifprofiles_short={}
for key in selectedMotSca_short:
	mm=key.split()[1]
	ss=int(key.split()[2])
	smoothed, norm= cleverMachine_petr(motifsshort[mm],property_table,AAs)
	motifprofiles_short[mm+" "+str(ss)]=smoothed[ss]


print "load motifs metapredictor"
motifslong={}
motifslong=create_dic_fromFasta("./data/motifsANDncMotifs30cdhit.txt")
motifprofiles_long={}
for key in selectedMotSca_long:
	mm=key.split()[1]
	ss=int(key.split()[2])
	smoothed, norm= cleverMachine_petr(motifslong[mm],property_table,AAs)
	motifprofiles_long[mm+" "+str(ss)]=smoothed[ss]
	
	

print "calculate Features"
super_model = joblib.load("./data/SuperModel/FinalRBFModel.model")
subprocess.call("python FilterMulti.py",shell=True) #creates feature vector for prediction

featdic={}
infile=open("./outputs/FinalFeatures.txt","r")
for line in infile:

	featvec=numpy.zeros([1, 60])
	for i in range(2,len(line.split())):
		featvec[0,i-2]=line.split()[i]
		#In [168]: pred_model.classes_
		#Out[168]: array([-1,  1])
		featdic[line.split()[0]]=super_model.predict_proba(featvec)[0][1]
		#featdic[line.split()[0]]=0.56

print featdic

fastafile="./data/inseq.fasta"
print "read fastafile ",fastafile
infile=open(fastafile,"r")
posRBP={}
name=""
for line in infile:
	if line[0]==">":
		name=line.split()[0]
	else:
		posRBP[name]=line
		name=""

#test=open("./outputs/FinalFeatures.txt","r").readline()


HMMERpresent=0
multio=open("./outputs/catRAPID_s_predictions.txt","w")
for RBP in posRBP:
	print RBP
	multio.write("--------------\n\n")	
	multio.write(RBP+"\n\n")
	print "check HMMER"
	o=open("./outputs/HMMER.fasta","w")
	o.write(">"+str(RBP)+"\n")
	o.write(posRBP[RBP]+"\n")
	o.close()
	subprocess.call("hmmscan --cut_ga --noali --domtblout ./outputs/fragmentation.log /mnt/fast/database/Pfam_hmm/Pfam-A.hmm ./outputs/HMMER.fasta > ./outputs/hmm_log.txt",shell=True)

	#subprocess.call("hmmscan --cut_ga --noali --domtblout ./outputs/fragmentation.log /mnt/large/databases/Pfam_hmm/Pfam-A.hmm ./outputs/HMMER.fasta > ./outputs/hmm_log.txt",shell=True)

	infile=open("./outputs/fragmentation.log","r")
	small={}
	for line in infile:
		if line[0]!="#":
			pfamid=line.split()[1].split(".")[0]
			small[pfamid+" "+line.split()[17]+" "+line.split()[18]]=line.split()[17]+" "+line.split()[18]+" "+line.split()[0] #id, position from to, domname
	infile.close()
	print small
	o=open("./outputs/HMMER.results","w")
	final={}
	if small!=[]: 
		for a in small:
			if (a.split()[0] in cpfam) or(a.split()[0] in ncpfam):#or(a.split()[0] in putpfam):
				o.write(a.split()[0]+" "+a.split()[1]+" "+a.split()[2]+" "+small[a].split()[2]+"\n")
				HMMERpresent=1
				final[a]=small[a]
			if (a.split()[0] in addpfam):
				HMMERpresent=1

	o.close()
	if os.stat("./outputs/HMMER.results").st_size == 0:
		o=open("./outputs/HMMER.results","w")
		o.write("Nodomain")
		o.close()

	check=[]
	
	if (featdic[RBP]>=0.5) or (HMMERpresent== 1):

		result=[]
		sequence=posRBP[RBP]
		RBPlength=len(sequence)
		posRBP_profiles={}
		for key in posRBP:
			smoothed, norm= cleverMachine_petr(posRBP[key],property_table,AAs)
			posRBP_profiles[key]=smoothed
		#load models for classical, nonclassical and putative prediction
		#print "load models for classical, nonclassical and putative prediction"

		#C_model = joblib.load("./data/Cmodel/FinalModel.model")	
		#scaletta={}
		#infile=open("./data/scaletta_meta/C_scaletta_0.75.txt","r")#49 0.4 ENSP00000236228_340-419 74 0.65511770892561
		#for line in infile: 
		#	if line[0]!="#":
		#		motif=line.split()[2]
		#		scale=line.split()[3]
		#		scaletta[motif+" "+scale]=line.split()[4]+" "+line.split()[1] #[motif scale]=cov corr
		#totaldata = numpy.zeros([1, len(scaletta)])		
		#totaldata= features(RBP,RBPlength,posRBP_profiles,scaletta,motifprofiles_long)
		#predC=C_model.predict_proba(totaldata)[0][1]
		#o.write("ClassicalScore "+str(round(pred,2))+"\n")
		

		#NC_model = joblib.load("./data/NCmodel/FinalModel.model")
		#scaletta={}
		#infile=open("./data/scaletta_meta/NC_scaletta_0.75.txt","r")#49 0.4 ENSP00000236228_340-419 74 0.65511770892561
		#for line in infile: 
		#	if line[0]!="#":
		#		motif=line.split()[2]
		#		scale=line.split()[3]
		#		scaletta[motif+" "+scale]=line.split()[4]+" "+line.split()[1] #[motif scale]=cov corr
		#totaldata = numpy.zeros([1, len(scaletta)])		
		#totaldata= features(RBP,RBPlength,posRBP_profiles,scaletta,motifprofiles_long)
		#predNC=NC_model.predict_proba(totaldata)[0][1]
		#o.write("NClassicalScore "+str(round(pred,2))+"\n")


		#P_model = joblib.load("./data/Pmodel/FinalModel.model")
		#scaletta={}
		#infile=open("./data/scaletta_meta/P_scaletta_0.75.txt","r")#49 0.4 ENSP00000236228_340-419 74 0.65511770892561
		#for line in infile: 
		#	if line[0]!="#":
		#		motif=line.split()[2]
		#		scale=line.split()[3]
		#		scaletta[motif+" "+scale]=line.split()[4]+" "+line.split()[1] #[motif scale]=cov corr
		#totaldata = numpy.zeros([1, len(scaletta)])		
		#totaldata= features(RBP,RBPlength,posRBP_profiles,scaletta,motifprofiles_long)
		#predP=P_model.predict_proba(totaldata)[0][1]
		#o.write("PlassicalScore "+str(round(pred,2))+"\n")

		
		#if featdic[RBP]>=0.5:	
		#	multio.write("Overall Interaction Score is "+str(round(featdic[RBP],2))+"\n\n")
		#else:
		#	multio.write("only HMMER region detected\n\n")


		print "load scan-model"
		pred_model=joblib.load("./data/ScanModel/FinalRBFscanModel.model")
	

		print "create prediction profile....."
		profile=open("./outputs/profile.txt","w")

		for centro in range(0,RBPlength):
			featVec=[]

			#realStart=int(a.split()[1])-3 #to adapt to smoothing profile
			#realEnd=int(a.split()[2])+3#to adapt to smoothing profile
			for couple in selectedMotSca_short:
				
				RRM=couple.split()[1]
				scale=int(couple.split()[2])
				currentRRMprof=motifprofiles_short[RRM+" "+str(scale)]
				RRMlength=len(motifsshort[RRM])#len(currentRRMprof)
				THEcorrelation=-1.0
				if (RRMlength+6 <= RBPlength-6): ###check if domain fits RBP
					gogo=0
					##centering RRM
					if (RRMlength% 2)== 0: ##even
						start=int(centro-int(RRMlength/2))
						end=int(centro+int(RRMlength/2))
					else:
						start=int(centro-(RRMlength/2))
						end=int(centro+(RRMlength/2)+1)

					#print start,end, end-start,RRMlength

					posseq=""
					if (start>=0)and(end<=RBPlength):

						if (start>=3) and (end+3<=RBPlength):
							posseq=posRBP[RBP][start-6:end+6]
						elif (start==2) and (end+4<=RBPlength):
							posseq=posRBP[RBP][start-5:end+7]
						elif (start==1) and (end+5<=RBPlength):
							posseq=posRBP[RBP][start-4:end+8]
						elif (end>=RBPlength) and (start>3):
							posseq=posRBP[RBP][start-12:end]
						elif (start==0) and (end+6)<=RBPlength:
							posseq=posRBP[RBP][start:end+12]
						elif (start>=3) and (end+2<=RBPlength):								
							posseq=posRBP[RBP][start-7:end+5]
						elif (start>=3) and (end+1<=RBPlength):								
							posseq=posRBP[RBP][start-8:end+4]
						elif (start>=3) and (end==RBPlength):								
							posseq=posRBP[RBP][start-12:end]			
						else:
							print "WTF"
						THEcorrelation=-1.0
						if posseq!="":
							THEcorrelation=-1.0
							smoothed, norm= cleverMachine_petr(posseq,property_table,AAs)
							posseq_profile=smoothed[scale]
							#posset.write(RBP+" "+couple+" "+posseq+"\n")
							THEcorrelations= numpy.zeros([1,6])


							if (len(posseq_profile)%2)==0: ##RBP even
								if (len(currentRRMprof)%2)==0: ##RD even
									THEstartpoint=(len(posseq_profile)/2)-(len(currentRRMprof)/2)
								else:
									THEstartpoint=(len(posseq_profile)/2)-(int(len(currentRRMprof)/2 +0.5))
							else:
								if (len(currentRRMprof)%2)==0: ##RD even
									THEstartpoint=(int(len(posseq_profile)/2+0.5))-(len(currentRRMprof)/2)
								else:
									THEstartpoint=(int(len(posseq_profile)/2+0.5))- (int(len(currentRRMprof)/2+0.5))


							if (THEstartpoint>=3):
								for sw in range(-3,3):							
									THEcorrelations[0,sw+3]=numpy.corrcoef(currentRRMprof,posseq_profile[THEstartpoint+sw:THEstartpoint+len(currentRRMprof)+sw])[0, 1]
								THEcorrelation=THEcorrelations.mean() 
							elif (THEstartpoint<3) and (THEstartpoint>0):
								if (THEstartpoint+len(currentRRMprof)+3)<=len(posseq_profile):
									if THEstartpoint==2:
										for sw in range(-2,4):						
											THEcorrelations[0,sw+3]=numpy.corrcoef(currentRRMprof,posseq_profile[THEstartpoint+sw:THEstartpoint+len(currentRRMprof)+sw])[0, 1]
									if THEstartpoint==1:
										for sw in range(-1,5):						
											THEcorrelations[0,sw+3]=numpy.corrcoef(currentRRMprof,posseq_profile[THEstartpoint+sw:THEstartpoint+len(currentRRMprof)+sw])[0, 1]
									THEcorrelation=THEcorrelations.mean() 
								else:
									THEcorrelations= numpy.zeros([1,5])
									for sw in range(-2,2):						
										THEcorrelations[0,sw+3]=numpy.corrcoef(currentRRMprof,posseq_profile[THEstartpoint+sw:THEstartpoint+len(currentRRMprof)+sw])[0, 1]
								THEcorrelation=THEcorrelations.mean() 
							else:
								THEcorrelation=-1.0
						else:
							THEcorrelation=-1.0
					else:
						THEcorrelation=-1.0
				else:
					THEcorrelation=-1.0
				featVec.append(round(THEcorrelation,2))

			#print centro,featVec

		
			count=0
			output="1 "
			for item in featVec:
				output+=str(item)+" "

			if output!="1 " and len(featVec)==60:
				print >>testseq, output

			lab=pred_model.predict_proba(featVec)[0]
			#if lab[1]>0.5:
			#	print "YES"
			result.append(lab[1])
			profile.write(str(centro)+" "+str(lab[0])+" "+str(lab[1])+"\n")	

		profile.close()
		print "smoothing"		
		for repe in range(10):
			result.insert(0,0.19179367326621416)
		for repe in range(10):
			result.append(0.19179367326621416)


		oo=open("./outputs/smoothy21.txt","w")
		for i in range(10,len(result)-10-1):
			sub=result[i-10:i+11]
			sum(sub)
			oo.write(str(i-10)+" "+str(sum(sub)/len(sub))+"\n")
		oo.close()



		#print "profile plotting..."
		#subprocess.call("cat profilePlotter2.R | R --slave --vanilla ",shell=True)	
		
		
                #if featdic[RBP]>=0.5:  
                #       multio.write("Overall Interaction Score is "+str(round(featdic[RBP],2))+"\n\n")
                #else:
                #       multio.write("only HMMER region detected\n\n")

		if featdic[RBP]>=0.5:		
			multio.write("Overall Interaction Score is "+str(round(featdic[RBP],2))+"\n\n")
			print "Overall Interaction Score is "+str(round(featdic[RBP],2))
			if HMMERpresent==1:
				print "HMMER region present"
				subprocess.call("python getBindingZones_mischmasch.py "+sequence,shell=True)
			else:
				print "No HMMER region present"
				subprocess.call("python getBindingZones.py "+sequence,shell=True)
			if os.stat("./outputs/list.txt").st_size !=0:
				print "write detected binding regions..."
				multio.write("Detected binding regions \n")
				multio.close()
				subprocess.call("cat ./outputs/list.txt >> ./outputs/catRAPID_s_predictions.txt",shell=True)
				multio=open("./outputs/catRAPID_s_predictions.txt","a")
			else:
				print "no binding regions found..."
	                        multio.write("No individual binding regions found \n")
                                #multio.close()
                                #subprocess.call("cat ./outputs/list.txt >> ./outputs/catRAPID_s_predictions.txt",shell=True)
                                #multio=open("./outputs/catRAPID_s_predictions.txt","a")

		else:
			if HMMERpresent==1:
				print "score below 0.5 but HMMER present"
				multio.write("Overall Interaction Score is "+str(round(featdic[RBP],2))+"\n\n")
				o=open("./outputs/list.txt","w")
				#multio.write("only HMMER region detected\n\n")
				print final
				if final!={}:
					multio.write("Detected HMMER regions \n")
					
					for i in final:			
						#multio.write(">"+i.split()[0]+"_"+i.split()[1]+"_"+i.split()[2]+"\n")
						#multio.write(sequence[int(i.split()[1])-1:int(i.split()[2])]+"\n")
	                        	        print i
		                                dist= int(int(i.split()[2])-int(i.split()[1]))
        		                        print dist
	                	                print len(sequence)
	                        	        if dist>=51:
        	                        	        print " ok"
                	                        	o.write(">"+i.split()[0]+"_"+i.split()[1]+"_"+i.split()[2]+"\n")
	                	                        o.write(sequence[int(i.split()[1])-1:int(i.split()[2])]+"\n")
        	                	        else:
	                	                        print "dist not ok, allungare"
                	        	                mezz= int(dist/2)+int(i.split()[1])
		                                        print mezz
        	        	                        ((mezz-26)>=0)
                	        	                ((mezz+25)<=len(sequence))
                        	        	        if ( ((mezz-26)>=0) and ((mezz+25)<=len(sequence)) ):
		                                                print "bothsites ok"
        		                                        o.write(">"+i.split()[0]+"_"+i.split()[1]+"_"+i.split()[2]+"_extendedTo_"+str(mezz-26)+"_"+str(mezz+25)+"\n")
                		                                o.write(sequence[mezz-26:mezz+25]+"\n")
                        		                else:
	                        	                        if (mezz-26<=0):
        	                        	                        o.write(">"+i.split()[0]+"_"+i.split()[1]+"_"+i.split()[2]+"_extendedTo_"+str(0)+"_"+str(mezz+25+(mezz-26))+"\n")
                	                        	                o.write(sequence[0:mezz+25+(mezz-26)]+"\n")
                        	                        	elif (mezz+25>=len(sequence)):
	                                	                        o.write(">"+i.split()[0]+"_"+i.split()[1]+"_"+i.split()[2]+"_extendedTo_"+str(mezz-25-(mezz-26))+"_"+str(len(sequence))+"\n")
        	                                	                o.write(sequence[mezz-25-(mezz-26):len(sequence)]+"\n")
					o.close()
					multio.close()
        		                subprocess.call("cat ./outputs/list.txt >> ./outputs/catRAPID_s_predictions.txt",shell=True)
                		        multio=open("./outputs/catRAPID_s_predictions.txt","a")

			else:
				print "score below 0.5 and no HMMER present"
				multio.write("No individual binding regions found  \n")


	else:
		print "No binding"
		multio.write("Overall Interaction Score is "+str(round(featdic[RBP],2))+"\n")
		#multio.write("\nno binding regions detected\n")
	
	multio.write("\n")
	subprocess.call("rm ./outputs/list.txt",shell=True)
	subprocess.call("rm ./outputs/HMMER*",shell=True)
	subprocess.call("rm ./outputs/hmm_log.txt",shell=True)
	subprocess.call("rm ./outputs/profile.txt",shell=True)
	subprocess.call("rm ./outputs/smoothy21.txt",shell=True)


multio.close()



randomnumber=sys.argv[1]
recipient=sys.argv[2]
#randomnumber=20
recipient="carmenmaria.livi@crg.es"
print "send email to ", recipient
#subprocess.call("python SendEmail.py "+recipient+" "+str(randomnumber),shell=True)

import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import sys


randomnumber=sys.argv[1]
email=sys.argv[2]

print randomnumber
print email.split()[0]


filename = "./outputs/catRAPID_s_predictions.txt"
f = file(filename)
attachment = MIMEText(f.read())
attachment.add_header('Content-Disposition', 'attachment', filename=filename)           



smtpObj = smtplib.SMTP('smtp.gmail.com', 587)
smtpObj.ehlo()
smtpObj.starttls()
smtpObj.ehlo()
smtpObj.login("catrapid.crg@gmail.com","gianga19761976")


msg = MIMEMultipart('alternative')
msg.attach(attachment)
toEmail, fromEmail =email.split()[0], "catrapid.crg@gmail.com"
msg['Subject'] = 'catRAPID signature submission no.: '+str(randomnumber)
msg['From'] = fromEmail
body = 'Please find the predicted RNA-binding regions of your submitted proteins attached.'
content = MIMEText(body, 'plain')
msg.attach(content)

try:
	f=open(filename,"r")
	f.close()
	smtpObj.sendmail(fromEmail, toEmail, msg.as_string())
except:
	print "Error, no email has been send"






