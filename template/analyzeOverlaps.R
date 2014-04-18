data=read.table("outputs/knownoverlaps.txt")

#RBPid,"overlap",RRMid,scale, mystart,myend, corr,err,hRRMid, hstart,hstop
####extracts overlap statistics
for (RBPid in unique(data$V1))
{
	#plusscore=0
	#minusscore=0
	oshift3=oshift4=0
	overlapWithHmmer=0
	set=data[which(RBPid==data$V1),]
	omaxmean=NULL
	overlap=set[which(set$V2=="overlap"),]
	for (RRMid in unique(overlap$V3))
	{
		d1=overlap[which(RRMid==overlap$V3),]
		for (scale in unique(d1$V4))
		{
			d2=d1[which(scale==d1$V4),] ##on overlapping profiletyp and scale
			if (length(unique(d2$V5))==3){oshift3=oshift3+1} ##couns if overlap is shift3
			if (length(unique(d2$V5))>3){oshift4=oshift4+1} ##or if overlap is shift4
			omaxmean=rbind(omaxmean,max(d2$V7))
			overlapWithHmmer=length(unique(d2$V9)) ##how many hmmer motifs are matched by the subset
			
		}
		
	}
	
	
	noshift3=noshift4=0
	nooverlap=set[which(set$V2=="NO-overlap"),]
	nomaxmean=NULL
	for (RRMid in unique(nooverlap$V3))
	{
		nd1=nooverlap[which(RRMid==nooverlap$V3),]
		for (scale in unique(nd1$V4))
		{
			nd2=nd1[which(scale==nd1$V4),]
			if (length(unique(d2$V5))==3){noshift3=noshift3+1}
			if (length(unique(d2$V5))>3){noshift4=noshift4+1}
			nomaxmean=rbind(nomaxmean,max(nd2$V7))
			overlapWithHmmer=length(unique(nd2$V9))
		}
	}
	write(paste(RBPid,"|",oshift3,oshift4,mean(omaxmean),"|",noshift3,noshift4,mean(nomaxmean),sep=" "),file="outputs/knownstatistic.txt",append=TRUE)
}


######NOT clear#####################

data=read.table("outputs/knownstatistic.txt")
p1 <- hist(data$V3)                     # centered at 4
p2 <- hist(data$V4)                     # centered at 6
p3 <- hist(data$V7)                     # centered at 4
p4 <- hist(data$V8)  
plot( p1$mids,p1$counts/length(data$V3), col=1,type="l",xlab="number of overlaps",ylab="normalized frequency")  # first histogram
lines( p2$mids,p2$counts/length(data$V4), col=2) 
lines( p3$mids,p3$counts/length(data$V7), col=3,type="l" )  # first histogram
lines( p4$mids,p4$counts/length(data$V8), col=4) 
legend('topright',c('overlaps_shift3','overlaps_shift4','Noverlaps_shift3','Noverlaps_shift4'),fill = c(1,2,3,4), bty = 'n',)

p1=hist(data$V5)
p2=hist(data$V9)
plot( p2, col=rgb(0,0,1,1/4),xlab="maxmean")  # first histogram
plot( p1, col=rgb(1,0,0,1/4), add=T) 
legend('topright',c('corr_overlap','corr_noOverlap'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n',)


####analyse TP rate (First try...not ok)
data=read.table("outputs/knownstatistic.txt")
TP=TN=0
for (i in seq(max(na.omit(round(data$V5,2))), min(na.omit(round(data$V5,2))),-0.01))
{
	oset=data[which(round(data$V5,2)==i),]
	for (ii in seq(1,nrow(oset)))
	{
		#print(paste(i,ii,sep=" "))
		d1=oset[ii,]
		if (nrow(d1)>0)
		{
			if (!(is.na(d1$V5)))
			{
			
				if (d1$V5>=d1$V9)
				{
					TP=TP+1	
				}else{
					TN=TN+1
				}
			}
		}	
	}
	print(paste(i,TP,TN,sep=" "))
}


####################vary correlation and calculate cases above correlation-threshold for overlapping and non-overlapping events
TP=TN=0
overlap=data$V5
noverlap=data$V9

#result=rbind("#corrthres","num_overlap","num_nonOverlap")
result=NA
for (i in seq(min(na.omit(round(data$V5,2))), max(na.omit(round(data$V5,2))),0.01)) #variing corr-threshold
{
	print(paste("threshold ",i,length(which(overlap>=i)),length(which(noverlap>=i))))
	result=rbind(result,c(i,length(which(overlap>=i)),length(which(noverlap>=i)),round(length(which(overlap>=i))/length(which(noverlap>=i)),2))) ##counting overlapping and nonoverlapping correlations above
}

print(result)
> result
       [,1] [,2] [,3] [,4]
result   NA   NA   NA   NA
       0.74  124  144 0.86
       0.75  120  144 0.83
       0.76  118  142 0.83
       0.77  113  138 0.82
       0.78  103  131 0.79
       0.79  100  125 0.80
       0.80   85  107 0.79
       0.81   75   93 0.81
       0.82   66   64 1.03
       0.83   52   47 1.11
       0.84   42   32 1.31
       0.85   33   25 1.32
       0.86   29   18 1.61
       0.87   23   15 1.53
       0.88   17    9 1.89
       0.89   15    7 2.14
       0.90   11    5 2.20
       0.91    8    2 4.00
       0.92    7    2 3.50
       0.93    6    2 3.00
       0.94    5    2 2.50
       0.95    5    2 2.50
       0.96    4    2 2.00
       0.97    4    2 2.00
       0.98    4    2 2.00
       0.99    4    2 2.00
       1.00    3    2 1.50

####################vary error and calculate cases above error-threshold for overlapping and non-overlapping events
TP=TN=0
overlap=data$V5
noverlap=data$V9

#result=rbind("#corrthres","num_overlap","num_nonOverlap")
result=NA
for (i in seq(min(na.omit(round(data$V,2))), max(na.omit(round(data$V,2))),0.01)) #variing corr-threshold
{
	print(paste("threshold ",i,length(which(overlap>=i)),length(which(noverlap>=i))))
	result=rbind(result,c(i,length(which(overlap>=i)),length(which(noverlap>=i)),round(length(which(overlap>=i))/length(which(noverlap>=i)),2))) ##counting overlapping and nonoverlapping correlations above
}

print(result)






####################vary correlation and take top 1,2,3.... and calculate TP/FP rate and F

data=read.table("outputs/knownoverlaps.txt")
sorted=data[order(data$V7,decreasing = TRUE),]

#subset=data[which(data$V7>=0.8),]
#sorted2=subset[order(subset$V7,decreasing = TRUE),]
for (pick in seq(2,nrow(data),1))
#for (pick in seq(2,20,1))
{
	#print(pick)
	set=sorted[1:pick,]
	TP=length(which(set$V2=="overlap"))
	FP=length(which(set$V2=="NO-overlap"))
	print(paste(pick,TP,FP,TP/FP,sep=" "))
}


1] "2 1 1 1"
[1] "3 1 2 0.5"
[1] "4 1 3 0.333333333333333"
[1] "5 1 4 0.25"
[1] "6 1 5 0.2"
[1] "7 2 5 0.4"
[1] "8 2 6 0.333333333333333"
[1] "9 2 7 0.285714285714286"
[1] "10 3 7 0.428571428571429"
[1] "11 3 8 0.375"
[1] "12 3 9 0.333333333333333"
[1] "13 3 10 0.3"
[1] "14 3 11 0.272727272727273"
[1] "15 3 12 0.25"
[1] "16 3 13 0.230769230769231"
[1] "17 3 14 0.214285714285714"
[1] "18 3 15 0.2"
[1] "19 3 16 0.1875"
[1] "20 3 17 0.176470588235294"


############analyze overlaps per protein to check if TP protein or FP protein

data=read.table("outputs/overlaps.txt")
file.remove("test.txt")
write("#Corr RBPid OverlapnumOfMyDomain OverlapnumOfHmmerDomains NoOverlapNumOfMyDomain NoOverlapegNumOfHmmerDomains",file="test.txt",append=TRUE)
#RBPid,"overlap",RRMid,scale, mystart,myend, corr,err,hRRMid, hstart,hstop
####extracts overlap statistics
for (corr in seq(1,0.6,-0.01))
{
	print(corr)
	for (RBPid in unique(data$V1))
	{	
		set=data[(RBPid==data$V1),]
		ol=set[(which("overlap"==set$V2)),]
		nol=set[(which("NO-overlap"==set$V2)),]

	
		numOfMyDomain=0
		numOfHmmerDomains=0
		myProfiles=unique(ol$V3)
		for (profile in myProfiles)
		{	
			ss=unique(ol[which(ol$V3==profile & ol$V7>=corr),]$V4)
			numOfMyDomain=numOfMyDomain+length(ss)
		}
		HmmerProfiles=unique(ol$V9)
		for (Hprofile in HmmerProfiles)
		{	
			position=unique(ol[which(ol$V9==Hprofile),]$V10)
			numOfHmmerDomains=numOfHmmerDomains+length(position)
		}
	
		#print(paste(corr,numOfMyDomain,numOfHmmerDomains,sep=" "))

		NegNumOfMyDomain=0
		NegNumOfHmmerDomains=0
		nmyProfiles=unique(nol$V3)
		for (nprofile in nmyProfiles)
		{	
			nss=unique(nol[which(nol$V3==nprofile & nol$V7>=corr),]$V4)
			NegNumOfMyDomain=NegNumOfMyDomain+length(nss)
		}
		nHmmerProfiles=unique(nol$V9)
		for (nHprofile in nHmmerProfiles)
		{	
			nposition=unique(nol[which(nol$V9==nHprofile),]$V10)
			NegNumOfHmmerDomains=NegNumOfHmmerDomains+length(nposition)
		}
	
		#print(paste(corr,NegNumOfMyDomain,NegNumOfHmmerDomains,sep=" "))
		write(paste(corr,RBPid,numOfMyDomain,numOfHmmerDomains,NegNumOfMyDomain,NegNumOfHmmerDomains,sep=" "),file="test.txt",append=TRUE)
	}
}

##analyse
data=read.table("test.txt")

for (corr in seq(1,0.6,-0.01)){
#corr RBPid OverlapnumOfMyDomain OverlapnumOfHmmerDomains NoOverlapNumOfMyDomain NoOverlapegNumOfHmmerDomains
	corrdata=data[which(data$V1==corr),]
	TP=nrow(corrdata[which(corrdata$V3>=corrdata$V5),])
	FP=nrow(corrdata[which(corrdata$V3<corrdata$V5),])
	print(paste(corr,TP,FP,sep=" "))
}

[1] "1 133 23"
[1] "0.99 132 24"
[1] "0.98 132 24"
[1] "0.97 132 24"
[1] "0.96 132 24"
[1] "0.95 132 24"
[1] "0.94 131 25"
[1] "0.93 127 29"
[1] "0.92 126 30"
[1] "0.91 119 37"
[1] "0.9 108 48"
[1] "0.89 94 62"
[1] "0.88 85 71"
[1] "0.87 77 79"
[1] "0.86 67 89"
[1] "0.85 58 98"
[1] "0.84 52 104"
[1] "0.83 42 114"
[1] "0.82 37 119"
[1] "0.81 33 123"
[1] "0.8 29 127"
[1] "0.79 24 132"
[1] "0.78 25 131"
[1] "0.77 24 132"
[1] "0.76 24 132"
[1] "0.75 24 132"
[1] "0.74 22 134"
[1] "0.73 22 134"
[1] "0.72 20 136"
[1] "0.71 20 136"
[1] "0.7 20 136"
[1] "0.69 20 136"
[1] "0.68 20 136"
[1] "0.67 20 136"
[1] "0.66 20 136"
[1] "0.65 20 136"
[1] "0.64 20 136"
[1] "0.63 20 136"
[1] "0.62 20 136"
[1] "0.61 20 136"
[1] "0.6 20 136"







