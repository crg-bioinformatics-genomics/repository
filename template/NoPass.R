##creates file for each entry
####cat GetResultsstraightplots.R | R --slave --vanilla --args <RBPname> <RBPlen>
data=read.table("outputs/preparedFrag.log")

pfams=read.table("data/pfam.id")
RNApfams=unique(pfams$V3)
pfams=read.table("data/DNA.bind.fam.txt")
DNApfams=unique(pfams$V1)


protLen=as.integer(unique(data$V4))
print(protLen)
allRRM=unique(as.character(data$V1))
file.remove("outputs/noPass.txt")
jpeg(filename = "outputs/Hmmerprofiles.jpeg",   width = 1000, height = 400)#, units = "px", pointsize = 12,  quality = 75,   bg = "white", res = NA,  type = c("cairo", "Xlib", "quartz"))


plot(0,ylab="", xlab="amino acid",yaxt="n",xaxt="n",xlim=c(0,protLen),col="white",type="l",lwd = 2,ylim=c(0,length(allRRM)*0.1+0.5), frame.plot=FALSE)
legend("topleft",c("RNA-binding domain","DNA-binding domain","other"),text.col = c(10,4,8))

a=seq(0,protLen,by=10)
a[1]=1
axis(1,c(a))

#axis(1,c(a=seq(1,nrow(RBP),by=5)))
ypoint=0
##print RRMs aligned above the protein
for (RRMid in allRRM)
{
	DNAset=which(DNApfams==RRMid)
	RNAset=which(RNApfams==RRMid)
	if (length(DNAset)>0){
		write(RRMid,"outputs/noPass.txt",append=TRUE)
		ypoint=ypoint+0.05
		d1=data[which(data$V1==RRMid),]
		ranges=data.frame(cbind(d1$V5,d1$V6))
		for (sw in seq(1,nrow(ranges))){
			start=ranges[sw,1]
			stop=ranges[sw,2]
			print(paste(start,stop,sw,sep=" "))
			#xRRM=cbind(seq(max(sw),max(sw)+nrow(RRM)),ypoint)
			#lines(xRRM,col=color,lwd = 3)
			xRRM=cbind(seq(start,stop),ypoint)
			lines(xRRM,col=4,lwd = 3)
		}
		text((start+10),ypoint, RRMid, cex=0.7, pos=3, col=4)
	}
	if (length(RNAset)>0){
		write(RRMid,"outputs/noPass.txt",append=TRUE)
		ypoint=ypoint+0.05
		d1=data[which(data$V1==RRMid),]
		ranges=data.frame(cbind(d1$V5,d1$V6))
		for (sw in seq(1,nrow(ranges))){
			start=ranges[sw,1]
			stop=ranges[sw,2]
			print(paste(start,stop,sw,sep=" "))
			#xRRM=cbind(seq(max(sw),max(sw)+nrow(RRM)),ypoint)
			#lines(xRRM,col=color,lwd = 3)
			xRRM=cbind(seq(start,stop),ypoint)
			lines(xRRM,col=10,lwd = 3)
		}
		text((start+10),ypoint, RRMid, cex=0.7, pos=3, col=10)
	}
	if (length(RNAset)==0 & length(DNAset)==0){
		write(RRMid,"outputs/noPass.txt",append=TRUE)
		ypoint=ypoint+0.05
		d1=data[which(data$V1==RRMid),]
		ranges=data.frame(cbind(d1$V5,d1$V6))
		for (sw in seq(1,nrow(ranges))){
			start=ranges[sw,1]
			stop=ranges[sw,2]
			print(paste(start,stop,sw,sep=" "))
			#xRRM=cbind(seq(max(sw),max(sw)+nrow(RRM)),ypoint)
			#lines(xRRM,col=color,lwd = 3)
			xRRM=cbind(seq(start,stop),ypoint)
			lines(xRRM,col=8,lwd = 3)
		}
		text((start+10),ypoint, RRMid, cex=0.7, pos=3, col=8)
	}
}
dev.off()

