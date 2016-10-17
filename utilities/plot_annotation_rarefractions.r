#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)
outfile = args[1]
type = args[2]
remainder = args[seq(3,length(args))]
names = remainder[seq(1,length(remainder),2)]
errcol = remainder[seq(2,length(remainder),2)]

filex = substr(outfile,nchar(outfile)-2,nchar(outfile))
if(filex=="pdf") {
  pdf(outfile)
} else if(filex=="png") {
  png(outfile)
} else {
    stop("unsupported type of output file.")
}

#type = "gene"
#names = c('/Shared/Au/jason/Code/NEWFUZZ/tempall/data/gene_rarefraction.txt','/Shared/Au/jason/Code/NEWFUZZ/tempall/data/gene_full_rarefraction.txt')
#errcol=c("#FF000088","#0000FF88")
errwid=4
datawid=3
limitwid=4
limitcol="#00000044"

read_count = 0
gene_any_count_max = 0
for(name in names) {
  d1<-read.table(name)
  if(max(d1[,1])>read_count) { read_count = max(d1[,1]) }
  if(max(d1[,4])>gene_any_count_max) { gene_any_count_max = max(d1[,4]) }
}
par(mar=c(5,5,4,2))
plot(1,type="n",ylim=c(0,gene_any_count_max*1.1),xlim=c(0,read_count*1.1),ylab=paste("Median ",type," count",sep=""),xlab="Reads sampled",bty="n",cex.lab=1.5,cex.axis=1.5)
z= 0
for(name in names) {
  z = z+1
  d1<-read.table(name)
  gene_any_count_median = max(d1[,3])
  abline(h=gene_any_count_median,lwd=limitwid,lty=2,col=limitcol)
  lines(d1[,1],d1[,3],lwd=datawid,col=errcol[z])
  for(i in seq(1,length(d1[,1]))){
    segments(d1[i,1],d1[i,3],x1=d1[i,1],y1=d1[i,4],col=errcol[z],lwd=errwid)
    segments(d1[i,1],d1[i,3],x1=d1[i,1],y1=d1[i,2],col=errcol[z],lwd=errwid)
  }
  #points(cbind(d1[,1],d1[,4]),pch="-",cex=4,col=errcol[z])
  #points(cbind(d1[,1],d1[,2]),pch="-",cex=4,col=errcol[z])
}
dev.off()
