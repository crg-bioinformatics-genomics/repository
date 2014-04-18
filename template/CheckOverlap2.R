##creates file for each entry
####cat GetResultsstraightplots.R | R --slave --vanilla --args <RBPname> <RBPlen>

data=read.table("outputs/errorCorrRBP_RRM.txt")
RBPid=unique(as.character(data$V1))
allRRM=as.character(unique(data$V2))

file.remove("outputs/overlaps.txt")

write(paste("#RBPid RRMid PFAMid scale mystart myend maxcorr bestcorr minerr besterr",sep=""),file="outputs/overlaps.txt",append=TRUE)
for (RRMid in allRRM)
{
	d1=data[which(data$V2==RRMid),]
	ss=unique(d1$V4)
#	protx=seq(1,protLen)
	for (scale in ss)
	{
		set=d1[which(d1$V4==scale),]

		maxcorr=max(set$V6)
		bestcorr=unique(set$V7)
		minerr=min(set$V8)
		besterr=unique(set$V9)
		mystart=min(set$V5)
		myend=mystart+nrow(set)-1
		overlap=(nrow(set))
		pfamid=unique(set$V3)## perche' non cambia overlaps??
		write(paste(RBPid,RRMid,pfamid,scale, mystart,myend,maxcorr,bestcorr,minerr,besterr,sep=" "),file="outputs/overlaps.txt",append=TRUE)
	}
}






