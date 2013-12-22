#cat heatmap.R | R --slave --args score domain domanlength
library(gplots)
library(fields)
args <- commandArgs()
#print(args)
score=as.numeric(args[5])
#print(score)
domain=as.character(args[6])
#print(domain)
RRMlen=as.numeric(args[7])
#RRMlen=RRMlen/2

data=read.table("./outputs/RBP.out")
title="RBPprofile scan"

# HD image
png("./outputs/HD_image.png",width=1920,height=480,pointsize=20)
#png("./outputs/image.png",width=750,height=160)

#newdata=rbind(zerovector,data)
#bla=rbind(zerovector,halfdata)

if (score>0.62){
	subtitle=paste("score: ",score,"    significant binding region found (",domain,"*)",sep="")
	image.plot(as.matrix(t(data)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"blue","white","red"), zlim=c(-1, 1),horizontal=T,legend.width=1,legend.shrink=0.2,	
	smallplot=c(.7,.85,.3,.35),bigplot=c(.1,.9,.6,.7),main=title,sub=subtitle)
	image(as.matrix(t(data)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"blue","white","red"), zlim=c(-1, 1),ad=TRUE)
#	axis(1,at=seq(0,1,0.1),labels=seq(0,length(data),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)
	#axis(1,at=seq(0,1,0.1),labels=seq(1,length(data),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)
	axis(1,at=seq(0,1,0.1),labels=seq(as.integer(RRMlen/2),as.integer(length(data)+(RRMlen/2)),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)
#	axis(1,at=seq(0,1.1,0.1),labels=seq(1,as.integer(length(data)+(RRMlen/2)),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)

	}else{
	new=data*0.62
	subtitle=paste("low score:    NO significant binding region found (",domain,"*)",sep="")
	image.plot(as.matrix(t(new)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"blue","white","red"), zlim=c(-1, 1),horizontal=T,legend.width=1,legend.shrink=0.2,	
	smallplot=c(.7,.85,.3,.35),bigplot=c(.1,.9,.6,.7),main=title,sub=subtitle)
	image(as.matrix(t(new)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"blue","white","red"), zlim=c(-1, 1),ad=TRUE)
#	axis(1,at=seq(0,1,0.1),labels=seq(1,length(new),as.integer(length(new)/10)),lwd=0,lwd.ticks=1)
	axis(1,at=seq(0,1,0.1),labels=seq(as.integer(RRMlen/2),as.integer(length(data)+(RRMlen/2)),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)

}



dev.off()

# NORMAL image
# png("./outputs/image.png",width=1920,height=480,pointsize=20)
png("./outputs/image.png",width=900,height=240)

#newdata=rbind(zerovector,data)
#bla=rbind(zerovector,halfdata)

if (score>0.62){
	subtitle=paste("score: ",score,"    significant binding region found (",domain,"*)",sep="")
	image.plot(as.matrix(t(data)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"blue","white","red"), zlim=c(-1, 1),horizontal=T,legend.width=1,legend.shrink=0.2,	
	smallplot=c(.7,.85,.3,.35),bigplot=c(.1,.9,.6,.7),main=title,sub=subtitle)
	image(as.matrix(t(data)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"blue","white","red"), zlim=c(-1, 1),ad=TRUE)
#	axis(1,at=seq(0,1,0.1),labels=seq(0,length(data),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)
	#axis(1,at=seq(0,1,0.1),labels=seq(1,length(data),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)
	axis(1,at=seq(0,1,0.1),labels=seq(as.integer(RRMlen/2),as.integer(length(data)+(RRMlen/2)),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)
#	axis(1,at=seq(0,1.1,0.1),labels=seq(1,as.integer(length(data)+(RRMlen/2)),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)

	}else{
	new=data*0.62
	subtitle=paste("low score:    NO significant binding region found (",domain,"*)",sep="")
	image.plot(as.matrix(t(new)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"blue","white","red"), zlim=c(-1, 1),horizontal=T,legend.width=1,legend.shrink=0.2,	
	smallplot=c(.7,.85,.3,.35),bigplot=c(.1,.9,.6,.7),main=title,sub=subtitle)
	image(as.matrix(t(new)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"blue","white","red"), zlim=c(-1, 1),ad=TRUE)
#	axis(1,at=seq(0,1,0.1),labels=seq(1,length(new),as.integer(length(new)/10)),lwd=0,lwd.ticks=1)
	axis(1,at=seq(0,1,0.1),labels=seq(as.integer(RRMlen/2),as.integer(length(data)+(RRMlen/2)),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)

}



dev.off()

#cat heatmap.R | R --slave --args score domain
#library(gplots)
#library(fields)
#args <- commandArgs()
#print(args)
#score=as.numeric(args[5])
#print(score)
#domain=as.character(args[6])
#print(domain)

#data=read.table("./outputs/RBP.out")
#title="RBProfile scan"
#png("./outputs/image.png",width=1920,height=480,pointsize=20)

#if (score>0.62){
#	subtitle=paste("score: ",score,"    significant binding region found (*)",sep="")
#	image.plot(as.matrix(t(data)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"blue","white","red"), zlim=c(-1, 1),horizontal=T,legend.width=1,legend.shrink=0.2,	
#	smallplot=c(.7,.8,.3,.35),bigplot=c(.1,.9,.6,.7),main=title,sub=subtitle)
#	image(as.matrix(t(data)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"blue","white","red"), zlim=c(-1, 1),ad=TRUE)
#	axis(1,at=seq(0,1,0.1),labels=seq(0,length(data),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)

#	}else{
#	subtitle=paste("score: ",score,"    NO significant binding region found (*) ",sep="")
#	image.plot(as.matrix(t(data)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"lightblue","white","firebrick2"), zlim=c(-1, 1),horizontal=T,legend.width=1,legend.shrink=0.2,
#	smallplot=c(.7,.8,.3,.35),bigplot=c(.1,.9,.6,.7),main=title,sub=subtitle)
#	image(as.matrix(t(data)),frame.plot=F,yaxt='n',xaxt='n',col=colorpanel(10,"lightblue","white","firebrick2"), zlim=c(-1, 1),ad=TRUE)
#	axis(1,at=seq(0,1,0.1),labels=seq(0,length(data),as.integer(length(data)/10)),lwd=0,lwd.ticks=1)
#}



#dev.off()

