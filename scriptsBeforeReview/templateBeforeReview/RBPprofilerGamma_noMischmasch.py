from __future__ import division
import os
import sys

sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/local/lib/python2.7/site-packages/pip-1.1-py2.7.egg')
sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/site-packages/pip-1.1-py2.7.egg')
sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7')
qsys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/plat-linux2')
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
import matplotlib
import collections
import sys
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



#sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/local/lib/python2.7/site-packages/pip-1.1-py2.7.egg')
#sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/site-packages/pip-1.1-py2.7.egg')
#sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7')
#sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/plat-linux2')
#sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/lib-tk')
#sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/lib-old')
#sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/lib-dynload')
#sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/local/lib/python2.7/site-packages')
#sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/lib/python2.7/site-packages')
#sys.path.insert(1,'/mnt/fast/processor/algorithm-signature/signatureEnv/local/lib/python2.7/site-packages/IPython/extensions')

print sys.path


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
		rbpprof=RBP_profiles[RBP][scale]
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



fastafile="./data/inseq.fasta"
print "read fastafile ",fastafile

infile=open(fastafile,"r")
posRBP={}
name=""
for line in infile:
	if line!=">input_protein":
		posRBP["RBP"]=line



AAs = "ARNDCQEGHILKMFPSTWYVX"    
titles, property_table  = prof.parse_scales(  os.path.join("", "./scales/scales_generated.txt"), AAs )    
#print titles

print "load C_scales"
selectedMotSca={}
infilename = "./data/C_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca["C "+motif+" "+scale]=corr

print "load NC_scales"
infilename = "./data/NC_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca["NC "+motif+" "+scale]=corr

print "load P_scales"
infilename ="./data/P_scaletta_0.75.txt"
print infilename
ff=open(infilename,"r") ######16 0.3 ENSP00000304350_1208-1343 17 0.664549863466898
for line in ff:
	if line[0]!="#":
		#print line
		motif=line.split()[2]
		scale=line.split()[3]
		corr=line.split()[1]
		selectedMotSca["P "+motif+" "+scale]=corr

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



print "load motifs"
motifs={}
motifs=create_dic_fromFasta("./data/motifs.txt")
motifprofiles={}
for key in selectedMotSca:
	mm=key.split()[1]
	ss=int(key.split()[2])
	smoothed, norm= cleverMachine_petr(motifs[mm],property_table,AAs)
	motifprofiles[mm+" "+str(ss)]=smoothed[ss]
	

super_model = joblib.load("./data/SuperModel/FinalRFModel.model")
subprocess.call("python Filter.py",shell=True)
test=open("./outputs/FinalFeatures.txt","r").readline()
featvec=numpy.zeros([1, len(selectedMotSca)])
for i in range(1,len(test.split())):
	featvec[0,i-1]=test.split()[i]
#In [168]: pred_model.classes_
#Out[168]: array([-1,  1])
decision=super_model.predict_proba(featvec)[0][1]
print decision


HMMERpresent=0
for RBP in posRBP:

	print "check HMMER"
	o=open("./outputs/HMMER.fasta","w")
	o.write(">"+str(RBP)+"\n")
	o.write(posRBP[RBP]+"\n")
	o.close()
	#subprocess.call("hmmscan --cut_ga --noali --domtblout ./outputs/fragmentation.log /mnt/large/databases/Pfam_hmm/Pfam-A.hmm ./outputs/HMMER.fasta > ./outputs/hmm_log.txt",shell=True)
        subprocess.call("hmmscan --cut_ga --noali --domtblout ./outputs/fragmentation.log /mnt/fast/database/Pfam_hmm/Pfam-A.hmm ./outputs/HMMER.fasta > ./outputs/hmm_log.txt",shell=True)
        #/mnt/fast/database
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
	o.close()
	if os.stat("./outputs/HMMER.results").st_size == 0:
		o=open("./outputs/HMMER.results","w")
		o.write("Nodomain")
		o.close()

	check=[]
	tablefile=open("outputs/testtable.html","w")
	tablefile.write("<tbody>\n")
	for a in final: ##creates html table with: region-length pfamID scale corr
		#print key, rbpout[key]
		if a.split()[0] not in check:
			check.append(a.split()[0])
			tablefile.write("<tr>\n")	
			tablefile.write("<td> <h3><a href=\"http://pfam.xfam.org/family/"+str(a.split()[0])+"\">"+str(a.split()[0])+"</a></h3>  </td>\n") #pfamlink
			tablefile.write("<td>"+str(final[a].split()[2])+"</td></tr>\n") #domainname
	tablefile.write("</tbody>\n")
	tablefile.close()


	if (decision>=0.5) or (HMMERpresent== 1):
		#open file for prediction values
		print "open file for prediction values"
		o=open("./outputs/prediction.txt","w")
		if decision>=0.5:		
			o.write("#prediction score>0.5 suggests RNA-binding\n")
		else:
			o.write("#onlyHMMER\n")
		o.write("TotalpredictionScore "+str(round(decision,2))+"\n")

		result=[]
		sequence=posRBP[RBP]
		RBPlength=len(sequence)
		posRBP_profiles={}
		for key in posRBP:
			smoothed, norm= cleverMachine_petr(posRBP[key],property_table,AAs)
			posRBP_profiles[key]=smoothed
		#load models for classical, nonclassical and putative prediction
		print "load models for classical, nonclassical and putative prediction"

		C_model = joblib.load("./data/Cmodel/CFinalRFModel.model")	
		scaletta={}
		infile=open("./data/C_scaletta_0.75.txt","r")#49 0.4 ENSP00000236228_340-419 74 0.65511770892561
		for line in infile: 
			if line[0]!="#":
				motif=line.split()[2]
				scale=line.split()[3]
				scaletta[motif+" "+scale]=line.split()[4]+" "+line.split()[1] #[motif scale]=cov corr
		totaldata = numpy.zeros([1, len(scaletta)])		
		totaldata= features(RBP,RBPlength,posRBP_profiles,scaletta,motifprofiles)
		pred=C_model.predict_proba(totaldata)[0][1]
		o.write("ClassicalScore "+str(round(pred,2))+"\n")
		

		NC_model = joblib.load("./data/NCmodel/NCFinalRFModel.model")
		scaletta={}
		infile=open("./data/NC_scaletta_0.75.txt","r")#49 0.4 ENSP00000236228_340-419 74 0.65511770892561
		for line in infile: 
			if line[0]!="#":
				motif=line.split()[2]
				scale=line.split()[3]
				scaletta[motif+" "+scale]=line.split()[4]+" "+line.split()[1] #[motif scale]=cov corr
		totaldata = numpy.zeros([1, len(scaletta)])		
		totaldata= features(RBP,RBPlength,posRBP_profiles,scaletta,motifprofiles)
		pred=NC_model.predict_proba(totaldata)[0][1]
		o.write("NClassicalScore "+str(round(pred,2))+"\n")


		P_model = joblib.load("./data/Pmodel/PFinalRBFModel.model")
		scaletta={}
		infile=open("./data/P_scaletta_0.75.txt","r")#49 0.4 ENSP00000236228_340-419 74 0.65511770892561
		for line in infile: 
			if line[0]!="#":
				motif=line.split()[2]
				scale=line.split()[3]
				scaletta[motif+" "+scale]=line.split()[4]+" "+line.split()[1] #[motif scale]=cov corr
		totaldata = numpy.zeros([1, len(scaletta)])		
		totaldata= features(RBP,RBPlength,posRBP_profiles,scaletta,motifprofiles)
		pred=P_model.predict_proba(totaldata)[0][1]
		o.write("PlassicalScore "+str(round(pred,2))+"\n")
		o.close()

		

		print "load scan-model"
		pred_model=joblib.load("./data/ScanModel/FinalRBFscanModel.model")
	

		print "create prediction profile....."
		profile=open("./outputs/profile.txt","w")

		for centro in range(0,RBPlength):
			featVec=[]

			#realStart=int(a.split()[1])-3 #to adapt to smoothing profile
			#realEnd=int(a.split()[2])+3#to adapt to smoothing profile
			for couple in selectedMotSca:
				
				RRM=couple.split()[1]
				scale=int(couple.split()[2])
				currentRRMprof=motifprofiles[RRM+" "+str(scale)]
				RRMlength=len(motifs[RRM])#len(currentRRMprof)
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

		oo=open("./outputs/smoothy5.txt","w")
		for i in range(10,len(result)-10-1):
			sub=result[i-2:i+3]
			sum(sub)
			oo.write(str(i-10)+" "+str(sum(sub)/len(sub))+"\n")
		oo.close()
		oo=open("./outputs/smoothy11.txt","w")
		for i in range(10,len(result)-10-1):
			sub=result[i-5:i+6]
			sum(sub)
			oo.write(str(i-10)+" "+str(sum(sub)/len(sub))+"\n")
		oo.close()
		oo=open("./outputs/smoothy21.txt","w")
		for i in range(10,len(result)-10-1):
			sub=result[i-10:i+11]
			sum(sub)
			oo.write(str(i-10)+" "+str(sum(sub)/len(sub))+"\n")
		oo.close()

		print "profile plotting..."
		subprocess.call("cat profilePlotter2.R | R --slave --vanilla ",shell=True)	
		
		#http://s.tartaglialab.com/new_submission/catrapid_omics_protein/55053
		if decision>=0.5:		
			subprocess.call("python getBindingZones.py",shell=True)
		else:
			o=open("./outputs/list.txt","w")
			for i in final:			
				o.write(">"+i.split()[0]+"_"+i.split()[1]+"_"+i.split()[2]+"\n")
				o.write(sequence[int(i.split()[1])-1:int(i.split()[2])]+"\n")
			o.close()

	else:
		
		print "No binding"
		o=open("./outputs/prediction.txt","w")
		o.write("#prediction score<0.5 suggests no RNA-binding\n")
		o.write("TotalpredictionScore "+str(round(decision,2))+"\n")
		o.close()
		o=open("./outputs/list.txt","w")
		o.write(">"+RBP+"\n")
		o.write(posRBP[RBP]+"\n")
		o.close()

#pred_model.classes_
#Out[74]: array([-1,  1])


