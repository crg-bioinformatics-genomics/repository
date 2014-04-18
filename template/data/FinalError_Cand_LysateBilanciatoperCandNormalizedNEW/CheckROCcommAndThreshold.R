##cat CheckPvalue2.R | R --slave 
##check AUC between pos and neg, for each domain+scale+accum_correlation 
##
library(pROC)
domain=read.table("../domains.txt")
ddcount=nrow(domain)
print("read domains")
#for (ps in 10:20){
#sequ=seq(100,200,5)
#sequ=seq(-10,10)
sequ=seq(0,10)
file.remove("ROCaccumAndThresSubset0.txt")
write(paste("#corr domain scales auc threshold spec sens",sep=" "),file="ROCaccumAndThresSubset0.txt",append=T)
for (ps in sequ){
	corr=(ps/10)
	ppos=ps+14

	##corr=(ppos-100)/50
	print(paste("correlation",corr,sep=" "))
	print(paste("position in array",ppos,sep=" "))
	for (scales in 0:79)
	{
		print(scales)
		pname=paste("/mnt/large/carmen/Profile/Candidate_accum/Candidate",scales,".pos",sep="")
		pos=read.table(pname)
		nname=paste("/mnt/large/carmen/Profile/Lysate_subset0_accum/Lysate_subset0",scales,".neg",sep="")
#		nname=paste("/mnt/large/carmen/Profile/Lysate_subsetForCandidate0/Lysate_subsetForCandidate0",scales,".neg",sep="")
		neg=read.table(nname)

		for (dd in 1:ddcount) 
		{
			enter=domain$V1[dd]
			#check AUC between pos and neg
			
				TP=TN=FP=FN=0
			  ##################positives#####################
				#take all with the same domain
				#remove zero vectors
				PremoveZero=pos[which(pos$V25>pos$V26),]
				NremoveZero=neg[which(neg$V25>neg$V26),]
				#get the entry of the combination domain-scale for each protein
				l1=PremoveZero[which(PremoveZero$V2==as.character(enter)),]
				n1=NremoveZero[which(NremoveZero$V2==as.character(enter)),]
				#get accumcorrelation at correlation position ppos
				n2=n1[,ppos]
				l2=l1[,ppos]
				##make average of 
				#aven=(mean(n2)*length(n2)+mean(l2)*length(l2))/(length(n2)*length(l2))
				
				#l2=l1[which(l1$V3==scales),]
#				for (check in l2$V4)
				#l2=l1[,ppos]
				#TP=length(l2[which(l2>=aven)])
				#FN=length(l2[which(l2<aven)])
			 ##################negatives#####################
				#FP=length(n2[which(n2>=aven)])
				#TN=length(n2[which(n2<aven)])
				#count how many values
				plabels=seq(l2)
				if (length(n2)==1){nlabels=1}else{nlabels=seq(n2)}
				if (length(plabels)>0)
				{ if (length(nlabels)>0)
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
								write(paste(corr,enter,scales,auc,vv,sep=" "),file="ROCaccumAndThresSubset0.txt",append=T)
							}else{write(paste(corr,enter,scales,auc,"threshold infinite",sep=" "),file="ROCaccumAndThresSubset0.txt",append=T)}
						}else{write(paste(corr,enter,scales,auc,"threshold is NULL",sep=" "),file="ROCaccumAndThresSubset0.txt",append=T)}
	

						#write(paste(corr,enter,scales,auc,sep=" "),file="ROCaccumAndThres.txt",append=T)
					}else{write(paste(corr,enter,scales,"no negatives available",sep=" "),file="ROCaccumAndThresSubset0.txt",append=T)}				
				}else{write(paste(corr,enter,scales,"no positives available",sep=" "),file="ROCaccumAndThresSubset0.txt",append=T)}				
			}
		}
}

