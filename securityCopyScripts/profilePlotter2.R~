jpeg(paste("./outputs/plot.jpeg",sep=""),width = 900, height = 400, units = "px", pointsize = 12,quality = 90)
#par(mfrow = c(2,2)) 

HMMER=read.table("./outputs/HMMER.results")

smoo=read.table("./outputs/smoothy21.txt")
plot(smoo$V1,smoo$V2,type="l",xaxt='n',ylim=c(-0.1,1.15),xlab="Amino Acid Position",ylab="RNA-binding Propensity")

if (HMMER$V1!="Nodomain")
{
	shift=0.8
	color=1
	for (pfamID in unique(HMMER$V1))
	{
		set=HMMER[which(HMMER$V1==pfamID),]
		shift=shift+0.06
		color=color+1
		for (ii in seq(nrow(set)))
		{	
			#print(ii)
			line=set[ii,]
			#print(line)
			arrows(line$V2,shift, line$V3,shift, length = 0.25, angle = 0,lwd=2, col=color)
			

		}
		text(1,shift,label=line$V1, col=color)
		
	}
}
xa=seq(0,length(smoo$V1),10)
xa[1]=1
axis(1, at = xa, labels = TRUE, tick = TRUE)
arrows(0, 0.5, length(smoo$V1), 0.5, length = 0.25, angle = 0,col="blue")


rect(0, -1, 30, 1.2, col=rgb(1, 0, 0.1, 0.1),border=0)
rect( max(smoo$V1)-30, -1, max(smoo$V1), 1.2, col=rgb(1, 0, 0.1, 0.1),border=0)
dev.off()





jpeg(paste("./outputs/miniplot.jpeg",sep=""),width = 400, height = 300, units = "px", pointsize = 12)
data=read.table("outputs/prediction.txt")
orddata=data[order(data$V2, decreasing = T),]

label=NULL
color=NULL
for (i in orddata$V1)
{
	if (i=="PlassicalScore")
	{
		label=rbind(label,"Putative") 
		color=rbind(color,"cornflowerblue")
	}
	#if (i=="TotalpredictionScore")
	#{
	#	label=rbind(label,"Overall") 
	#	color=rbind(color,"firebrick1")
	#}
	if (i=="ClassicalScore")
	{
		label=rbind(label,"Classical") 
		color=rbind(color,"firebrick1")
	}
	if (i=="NClassicalScore")
	{
		label=rbind(label,"Non-classical") 
		color=rbind(color,"forestgreen")
	}
	if (i=="TotalpredictionScore")
	{
		totalScore=orddata[which(orddata$V1=="TotalpredictionScore"),]$V2
	}
}


#value=orddata$V2
value=orddata[which(orddata$V1!="TotalpredictionScore"),]$V2

#barplot(value, main="Propensity\nScores", ylim=c(0,1.15), col=color, ylab="Score", names.arg=label)
b=barplot(value, main=paste("Overall Interaction Score = ",totalScore,sep=" "), ylim=c(0,1.15), col=color, ylab="Interaction Propensity", names.arg=label)
text(x=b,y=value,labels=value,pos=3,col="black",cex=1.25)
dev.off()
