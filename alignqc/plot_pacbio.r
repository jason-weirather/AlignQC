#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)
if(length(args)<2) {
  stop("Must supply input and output\n",call.=FALSE)
}

# decide output type
filex = substr(args[2],nchar(args[2])-2,nchar(args[2]))
if(filex=="pdf") {
  pdf(args[2],bg="#FFFFFF",height=3.5,width=7)
} else if (filex=="png") {
  png(args[2],bg="#FFFFFF",width=480,height=240)
} else {
    stop("Unsupported type for output file.\n",call.=FALSE)
}

d<-read.table(args[1])
tot = length(d[,1])
biggest = max(d[,2])
par(mar=c(5,4.1,1,1))
plot(1,type="n",xlim=c(0,max(20,tot)+1),ylim=c(0,biggest),xlab="SMRT cell",ylab="Molecule count")
for(i in seq(1,tot,1)) {
  rect(-0.4+i,0,i+0.4,d[i,1],col="#777777")
  rect(-0.4+i,d[i,1],i+0.4,d[i,2],col="#FF0000")
}
dev.off()
