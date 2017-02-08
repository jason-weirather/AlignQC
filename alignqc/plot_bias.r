#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)
if(length(args)<2) {
  stop("Must supply input and output\n",call.=FALSE)
}

# decide output type
filex = substr(args[2],nchar(args[2])-2,nchar(args[2]))
if(filex=="pdf") {
  pdf(args[2],bg="#FFFFFF",width=7,height=4)
} else if (filex=="png") {
  png(args[2],bg="#FFFFFF",width=800,height=400)
} else {
    stop("Unsupported type for output file.\n",call.=FALSE)
}

filex = substr(args[1],nchar(args[1])-1,nchar(args[1]))
if(filex=='gz') {
  d<-read.table(gzfile(args[1]),row.names=1)
} else {
  d<-read.table(args[1],row.names=1)
}
myboxcol = '#FFDB58'
if(length(args)>2) {
  myboxcol = args[3]
}
par(mar=c(5,5,1,1))
axcex = 1.75
recwid = 3
boxplot(as.matrix(d),outline=FALSE,xlab="Position on Reference (5' to 3' %)",ylab="Coverage",lty=1,boxwex=1,staplecol="#77777788",whiskcol="#77777788",col=myboxcol,cex.axis=axcex,cex.lab=axcex,whisklwd=2,staplelwd=3,boxlwd=1,boxcol="#000000AA",names=seq(1,100,1),frame=FALSE,axes=FALSE)
axis(1,cex.axis=axcex,lwd=recwid)
axis(2,cex.axis=axcex,lwd=recwid)
dev.off()
