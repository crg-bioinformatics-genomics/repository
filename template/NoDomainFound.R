##creates file for each entry
####cat GetResultsstraightplots.R | R --slave --vanilla --args <RBPname> <RBPlen>
file.remove("outputs/noPass.txt")
write("NA",file="outputs/noDomainFound.txt")
jpeg(filename = "outputs/Hmmerprofiles.jpeg",   width = 1000, height = 400)#, units = "px", pointsize = 12,  quality = 75,   bg = "white", res = NA,  type = c("cairo", "Xlib", "quartz"))

plot(0,ylab="", xlab="amino acid",yaxt="n",main="NO domain found",xaxt="n",xlim=c(0,1),col="white",type="l",lwd = 2,ylim=c(0,1), frame.plot=FALSE)

dev.off()

