#cat heatmap.R | R --slave --args score domain
library(gplots)
library(fields)
args <- commandArgs()
print(args)
score=as.numeric(args[5])
print(score)
domain=as.character(args[6])
print(domain)

data=read.table("./outputs/RBP.out")
title="RBProfile scan"
png("./outputs/image.png",width=1920,height=480,pointsize=20)

if (score>0.62){
	subtitle=paste("score: ",score,"    significant binding region found (*)",sep="")
	image.plot(as.matrix(t(data)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"blue","white","red"), zlim=c(-1, 1),horizontal=T,legend.width=1,legend.shrink=0.2,	
	smallplot=c(.7,.8,.3,.35),bigplot=c(.1,.9,.6,.7),main=title,sub=subtitle)
	image(as.matrix(t(data)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"blue","white","red"), zlim=c(-1, 1),ad=TRUE)
	axis(1,at=seq(0,1,0.1),labels=seq(0,length(data),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)

	}else{
	subtitle=paste("score: ",score,"    NO significant binding region found (*) ",sep="")
	image.plot(as.matrix(t(data)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"lightblue","white","firebrick2"), zlim=c(-1, 1),horizontal=T,legend.width=1,legend.shrink=0.2,
	smallplot=c(.7,.8,.3,.35),bigplot=c(.1,.9,.6,.7),main=title,sub=subtitle)
	image(as.matrix(t(data)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"lightblue","white","firebrick2"), zlim=c(-1, 1),ad=TRUE)
	axis(1,at=seq(0,1,0.1),labels=seq(0,length(data),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)
}



dev.off()

