##creates file for each entry
####cat GetResultsstraightplots.R | R --slave --vanilla --args <RBPname> <RBPlen>
hmmer=read.table("outputs/preparedFrag.log")

pfams=read.table("data/pfam.id")
RNApfams=unique(pfams$V3)
pfams=read.table("data/DNA.bind.fam.txt")
DNApfams=unique(pfams$V1)

data=read.table("outputs/errorCorrRBP_RRM.txt")
RBPid=unique(as.character(data$V1))
allRRM=as.character(unique(data$V2))


protLen=as.integer(unique(hmmer$V4))
print(protLen)
RRMhmmer=unique(as.character(hmmer$V1))
file.remove("outputs/overlaps.txt")
#write("#RRMid scale mystart myend corr err hRRMid hstart hstop",file="outputs/overlaps.txt",append=TRUE)
#jpeg(filename = "outputs/overaps.jpeg",   width = 1000, height = 400)#, units = "px", pointsize = 12,  quality = 75,   bg = "white", res = NA,  type = c("cairo", "Xlib", "quartz"))


plot(0,ylab="", xlab="amino acid",yaxt="n",xaxt="n",xlim=c(0,protLen),col="white",type="l",lwd = 2,ylim=c(0,length(allRRM)*0.1), frame.plot=FALSE)
legend("topleft",c("RNA-binding domain","DNA-binding domain","other"),text.col = c(10,4,8))


a=seq(0,protLen,by=10)
a[1]=1
axis(1,c(a))

#axis(1,c(a=seq(1,nrow(RBP),by=5)))
ypoint=0

for (RRMid in allRRM)
{
	d1=data[which(data$V2==RRMid),]
	ss=unique(d1$V3)
	protx=seq(1,protLen)
	for (scale in ss)
	{

		set=d1[which(d1$V3==scale),]
		RBPname=paste("data/proteinProfiles/",RBPid,"_",scale,".txt",sep="")
		RBP=read.table(RBPname)
		RRMname=""
		RRMname=paste("data/domainProfiles/",RRMid,"_",scale,".txt",sep="")
		RRM=read.table(RRMname)
		for (line in set$V4)
		{
			corr=set[which(set$V4==line),]$V5
			err=set[which(set$V4==line),]$V7
			mystart=line
			myend=line+nrow(RRM)
			#print(paste(RRMid, scale,mystart,myend,sep=" "))		
#			write(paste("P ",RRMid, scale,mystart,myend,sep=" "),file="outputs/overlaps.txt",append=TRUE)
			for (hRRMid in RRMhmmer)
			{
				DNAset=which(DNApfams==hRRMid)
				RNAset=which(RNApfams==hRRMid)
				if (length(RNAset)>0)
				{
					#write(hRRMid,"outputs/overlaps.txt",append=TRUE)
					#ypoint=ypoint+0.05
					h=hmmer[which(hmmer$V1==hRRMid),]
					ranges=data.frame(cbind(h$V5,h$V6))
					for (hmmerrange in seq(1,nrow(ranges)))
					{
						hstart=ranges[hmmerrange,1]
						hstop=ranges[hmmerrange,2]
						#print(paste(hRRMid,hstart,hstop,hmmerrange,sep=" "))
						#write(paste("P ",RRMid, scale,mystart,myend,sep=" "),file="outputs/overlaps.txt",append=TRUE)
						#write(paste("H ",hRRMid, hstart,hstop,sep=" "),file="outputs/overlaps.txt",append=TRUE)

						print(paste(mystart,myend,sep=" "))
						if (((myend>=hstart)&(myend<=hstop)) | ((mystart>=hstart)&(mystart<=hstop)))
						{ 
							#print(paste("overlap",RRMid,hRRMid,sep=" "))
							write(paste(RBPid,"overlap",RRMid,scale, mystart,myend, corr,err,hRRMid, hstart,hstop,sep=" "),file="outputs/overlaps.txt",append=TRUE)
						}else{write(paste(RBPid,"NO-overlap",RRMid,scale, mystart,myend, corr,err,hRRMid, hstart,hstop,sep=" "),file="outputs/overlaps.txt",append=TRUE)}
						#xRRM=cbind(seq(max(sw),max(sw)+nrow(RRM)),ypoint)
						#lines(xRRM,col=color,lwd = 3)
						#xRRM=cbind(seq(start,stop),ypoint)
						#lines(xRRM,col=10,lwd = 3)
					}
					#text((start+10),ypoint, hRRMid, cex=0.7, pos=3, col=10)
				}else
				{
					#print("no RNA/DNA domain found")
					write(paste(RBPid,RRMid,"no RNA/DNA domain found by hmmer",sep=" "),file="outputs/noDomainFoundHmmer.txt",append=TRUE)
				}
			}
		}
	}
}






