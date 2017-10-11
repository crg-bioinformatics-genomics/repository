jpeg(paste("./outputs/plot.jpeg",sep=""),width = 1300, height = 500, units = "px", pointsize = 12,quality = 90)
par(mfrow = c(2,2)) 

HMMER=read.table("./outputs/HMMER.results")

data=read.table("./outputs/profile.txt")
plot(data$V1,data$V3,col="blue",type="l",xaxt='n',ylim=c(-0.1,1.1),xlab="amino acid",ylab="prediction score")
if (HMMER$V1!="Nodomain")
{
	arrows(HMMER$V2,1.0, HMMER$V3,1.0, length = 0.25, angle = 0,lwd=2, col="red")
}
axis(1, at = seq(1,length(data$V1),10), labels = TRUE, tick = TRUE)
arrows(0, 0.5, length(data$V1), 0.5, length = 0.25, angle = 0)


smoo=read.table("./outputs/smoothy5.txt")
plot(smoo$V1,smoo$V2,col="blue",type="l",xaxt='n',ylim=c(-0.1,1.1),xlab="amino acid",ylab="prediction score")
if (HMMER$V1!="Nodomain")
{
	arrows(HMMER$V2,1.0, HMMER$V3,1.0, length = 0.25, angle = 0,lwd=2, col="red")
}
axis(1, at = seq(1,length(smoo$V1),10), labels = TRUE, tick = TRUE)
arrows(0, 0.5, length(smoo$V1), 0.5, length = 0.25, angle = 0)


smoo=read.table("./outputs/smoothy11.txt")
plot(smoo$V1,smoo$V2,col="blue",type="l",xaxt='n',ylim=c(-0.1,1.1),xlab="amino acid",ylab="prediction score")
if (HMMER$V1!="Nodomain")
{
	arrows(HMMER$V2,1.0, HMMER$V3,1.0, length = 0.25, angle = 0,lwd=2, col="red")
}
axis(1, at = seq(1,length(smoo$V1),10), labels = TRUE, tick = TRUE)
arrows(0, 0.5, length(smoo$V1), 0.5, length = 0.25, angle = 0)

smoo=read.table("./outputs/smoothy21.txt")
plot(smoo$V1,smoo$V2,col="blue",type="l",xaxt='n',ylim=c(-0.1,1.1),xlab="amino acid",ylab="prediction score")
if (HMMER$V1!="Nodomain")
{
	arrows(HMMER$V2,1.0, HMMER$V3,1.0, length = 0.25, angle = 0,lwd=2, col="red")
}
axis(1, at = seq(1,length(smoo$V1),10), labels = TRUE, tick = TRUE)
arrows(0, 0.5, length(smoo$V1), 0.5, length = 0.25, angle = 0)
dev.off()


jpeg(paste("./outputs/miniplot.jpeg",sep=""),width = 400, height = 250, units = "px", pointsize = 12)
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
	if (i=="TotalpredictionScore")
	{
		label=rbind(label,"Overall") 
		color=rbind(color,"firebrick1")
	}
	if (i=="ClassicalScore")
	{
		label=rbind(label,"Classical") 
		color=rbind(color,"forestgreen")
	}
	if (i=="NClassicalScore")
	{
		label=rbind(label,"Non-classical") 
		color=rbind(color,"gold")
	}
}


value=orddata$V2

barplot(value, main="Propensity\nScores", ylim=c(0,1.15), col=color, ylab="Score", names.arg=label)
b=barplot(value, main="Propensity\nScores", ylim=c(0,1.15), col=color, ylab="Score", names.arg=label)
text(x=b,y=value,labels=value,pos=3,col="black",cex=1.25)
dev.off()
