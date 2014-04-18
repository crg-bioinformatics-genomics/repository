##cat CheckPvalue2.R | R --slave 
##check AUC between pos and neg, for each domain+scale+accum_correlation 
##
library(pROC)
domain=read.table("../domains.txt")

print("read domains")
file.remove("ROCminimumErrorSubset0.txt")
write(paste("#domain","scale","auc","bestThres","sens","spec",sep=" "),file="ROCminimumErrorSubset0.txt",append=T)
scalefile=read.table("/home/carmen/Documents/CRG/Profile_motifs/new2014/results/DomainsAndBestScale_pval0.1.NEW")
ddcount=nrow(scalefile)



for (item in 1:ddcount)
{	line=scalefile[item,]
	print(line)
	domain=line$V1

	scale=line$V2

#	pname=paste("/mnt/large/carmen/Profile/ErrorAUC/FinalError_Cand_LysateBil0/error",domain,scale,".pos",sep="")
#	print(pname)
#	pos=read.table(pname)
#	nname=paste("/mnt/large/carmen/Profile/ErrorAUC/FinalError_Cand_LysateBil0/error",domain,scale,".neg",sep="")
#	neg=read.table(nname)
#	print(nname)
	

	pname=paste("/mnt/large/carmen/Profile/ErrorAUC/FinalErrorNormalized_Cand_LysateSubsetForCandidates0NEW/error",domain,scale,".pos",sep="")
	print(pname)
	pos=read.table(pname)
	nname=paste("/mnt/large/carmen/Profile/ErrorAUC/FinalErrorNormalized_Cand_LysateSubsetForCandidates0NEW/error",domain,scale,".neg",sep="")
	neg=read.table(nname)
	print(nname)


	#coords(smooth.roc, x, input=c("specificity","sensitivity"), ret=c("specificity", "sensitivity"), as.list=FALSE,drop=TRUE, best.method=c("youden"))
	
	#check AUC between pos and neg	
	TP=TN=FP=FN=0
	l1=pos[which(pos$V2==as.character(domain)),]
	n1=neg[which(neg$V2==as.character(domain)),]
		#get accumcorrelation at correlation position ppos
	n2=n1[,4] #ENSP00000366843 ENSP00000318195_509-559 1 2.2
	l2=l1[,4]
	plabels=seq(l2)
	nlabels=seq(n2)
	if (length(plabels)>0)
	{ 
		if (length(nlabels)>0)
		{	
			for (i in plabels){plabels[i]="Good"}
			for (i in nlabels){nlabels[i]="Poor"}
			#posRocdata=cbind(plabels,l2)							
			test<-data.frame(label=c(plabels,nlabels),values=c(l2,n2))
			test2=test[order(test$values,decreasing=T),]
			#test=predictions[(order(as.numeric(predictions[,2]),decreasing=T)),]
			outcome=roc(test2$label,test2$value,plot=F,levels=c('Good','Poor'))#,direction=">")
			#outcome=roc(predictions[order(predictions[,2]),1],as.numeric(predictions[order(predictions[,2]),2]),plot=T,print.auc=F,percent=F)
			auc=outcome$auc
			values=coords(outcome, x="best",best.method="y",ret=c("threshold","specificity","sensitivity"),as.list=T)
			if (!is.null(values$threshold))
			{
				if (is.finite(values$threshold))
				{
					vv=paste(round(values$threshold,2),round(values$specificity,2),round(values$sensitivity,2),sep=" ")
					write(paste(domain,scale,auc,vv,sep=" "),file="ROCminimumErrorSubset0.txt",append=T)
				}else{write(paste(domain,scale,auc,"is infinite",sep=" "),file="ROCminimumErrorSubset0.txt",append=T)}
			}else{write(paste(domain,scale,auc,"is null",sep=" "),file="ROCminimumErrorSubset0.txt",append=T)}
		}else{write(paste(domain,scale,"no negatives available",sep=" "),file="ROCminimumErrorSubset0.txt",append=T)}				
		
	}else{write(paste(enter,scale,"no positives available",sep=" "),file="ROCminimumErrorSubset0.txt",append=T)}				
	
}
