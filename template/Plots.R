##creates file for each entry

args <- commandArgs()
print(args)

protLen=as.integer(args[5])
print(protLen)

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
	protx=seq(1,nrow(RBP))
	print(nrow(RBP))	
	jpeg(filename,   width = 1000, height = 400)
	plot(protx,RBP$V1,type="l",xlab="amino acid",ylab="property range")
	#plot(0,ylab="", xlab="amino acid",yaxt="n",xaxt="n",xlim=c(0,protLen),col="white",type="l",lwd = 2,ylim=c(0,length(allRRM)*0.1), frame.plot=FALSE)
	for (i in seq(d$V5,d$V6))
	{
		domainx=seq(i,i+nrow(RRM)-1)
		lines(domainx,RRM$V1,col=runvar)
	}
	dev.off()
}


jpeg(filename = "outputs/profiles.jpeg",   width = 1000, height = 400)#, units = "px", pointsize = 12,  quality = 75,   bg = "white", res = NA,  type = c("cairo", "Xlib", "quartz"))
plot(0,ylab="", xlab="amino acid",yaxt="n",xaxt="n",xlim=c(0,protLen),col="white",type="l",lwd = 2,ylim=c(0,length(allRRM)*0.1), frame.plot=FALSE)
a=seq(0,protLen,by=10)
a[1]=1
axis(1,c(a))

#axis(1,c(a=seq(1,nrow(RBP),by=5)))
color=1
ypoint=0
##print RRMs aligned above the protein
for (RRMid in allRRM)
{
	set=data[which(data$V2==RRMid),] 
	ss=unique(set$V4)
	color=color+1
	ypoint=ypoint+0.05
	
	RRMname=paste("data/domainProfiles/",RRMid,"_",ss,".txt",sep="")
	RRM=read.table(RRMname)

	xRRM=cbind(seq(set$V5,set$V6+nrow(RRM)),ypoint)
	lines(xRRM,col=color,lwd = 3)
		
	text((set$V6+(nrow(RRM)/2)),ypoint, set$V3, cex=0.7, pos=3, col=color)
}
dev.off()

