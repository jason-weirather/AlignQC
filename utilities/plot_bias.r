#!/opt/R/3.2.1/bin/Rscript

args=commandArgs(trailingOnly=TRUE)
if(length(args)<2) {
  stop("Must supply input and output\n",call.=FALSE)
}

# decide output type
filex = substr(args[2],nchar(args[2])-2,nchar(args[2]))
if(filex=="pdf") {
  pdf(args[2],bg="#FFFFFF",width=10,height=5)
} else if (filex=="png") {
  png(args[2],bg="#FFFFFF",width=800,height=400)
} else {
    stop("Unsupported type for output file.\n",call.=FALSE)
}

d<-read.table(gzfile(args[1]),row.names=1)
par(mar=c(5,5,1,1))
boxplot(t(as.matrix(d)),outline=FALSE,xlab="Position on covered portion of reference 5'-3' (%)",ylab="Fraction of transcripts covered",lty=1,boxwex=0.75,whiskcol="#77777733",col="#FFDB58",cex.axis=1.5,cex.lab=1.5,whisklwd=2,boxlwd=1,boxcol="#000000AA")
dev.off()
