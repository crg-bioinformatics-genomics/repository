##creates file for each entry

data=read.table("outputs/errorCorrRBP_RRM.Final")

RBPid=unique(as.character(data$V1))

allRRM=as.character(unique(data$V2))

runvar=1
for (RRMid in allRRM)
{
	d=data[which(data$V2==RRMid),]
	scale=unique(d$V4)
	runvar=runvar+1

	RBPname=paste("data/proteinProfiles/",RBPid,"_",scale,".txt",sep="")
	RBP=read.table(RBPname)
	RRMname=""
	RRMname=paste("data/domainProfiles/",RRMid,"_",scale,".txt",sep="")
	RRM=read.table(RRMname)
		
	filename=paste("outputs/",RRMid,"_",scale,".jpeg",sep="")
	print(filename)
	jpeg(filename,   width = 1000, height = 400)
	#	jpeg(filename)
	protx=seq(1,nrow(RBP))
	print(nrow(RBP))
	plot(protx,RBP$V1,type="l",xlab="amino acid",ylab="property range")
	#plot(0,ylab="", xlab="amino acid",yaxt="n",xaxt="n",xlim=c(0,protLen),col="white",type="l",lwd = 2,ylim=c(0,length(allRRM)*0.1), frame.plot=FALSE)
	for (i in seq(d$V5,d$V6))
	{
		domainx=seq(i,i+nrow(RRM)-1)
		lines(domainx,RRM$V1,col=runvar)
	}
		
	

	dev.off()
	
}
