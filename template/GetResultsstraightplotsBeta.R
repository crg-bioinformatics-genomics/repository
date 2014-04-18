##creates file for each entry
####cat GetResultsstraightplots.R | R --slave --vanilla --args <RBPname> <RBPlen>
data=read.table("outputs/errorCorrRBP_RRM.Final")
args <- commandArgs()
print(args)

#RBPid=args[5]
protLen=as.integer(args[5])
print(protLen)
allRRM=unique(as.character(data$V2))

jpeg(filename = "outputs/profiles.jpeg",   width = 1000, height = 400)#, units = "px", pointsize = 12,  quality = 75,   bg = "white", res = NA,  type = c("cairo", "Xlib", "quartz"))
#
#jpeg("outputs/profiles.jpeg", units="cm", width=10, height=5, res=100)

#RBPname=paste("proteinProfiles/",RBPid,"_",0,".txt",sep="") #get protein profile for scale 0 -->only to have protein length
#RBP=read.table(RBPname)
	
#protplot=cbind(seq(1,nrow(RBP)),0)
#plot(0,xlab="amino acid",yaxt="n",xaxt="n",xlim=c(0,nrow(RBP)),col=5,type="l",lwd = 2,ylim=c(0,length(allRRM)/10+0.1), frame.plot=FALSE)

##plot proteinlength
#protplot=cbind(seq(1,protLen),0)

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

