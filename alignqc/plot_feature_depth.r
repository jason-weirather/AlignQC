#!/usr/bin/env Rscript

# Read the inputs
# 1. genome_depth_bed
# 2. exon_depth_bed
# 3. intron_depth_bed
# 4. intergenic_depth_bed
# 5. output_image_file
args=commandArgs(trailingOnly=TRUE)
if(length(args)<5) {
  stop("Must supply inputs and output \n",call.=FALSE)
}
## decide output type
filex = substr(args[5],nchar(args[5])-2,nchar(args[5]))
if(filex=="pdf") {
  pdf(args[5],bg="#FFFFFF")
} else if (filex=="png") {
  png(args[5],bg="#FFFFFF")
} else {
    stop("Unsupported type for output file.\n",call.=FALSE)
}

dg<-read.table(gzfile(args[1]))
dexon<-read.table(gzfile(args[2]))
dintron<-read.table(gzfile(args[3]))
dinter<-read.table(gzfile(args[4]))

logtrans<-function(num) {
 if(num==1) {
   return(1)
 }
 return(log(num,2)+1)
}
untrans<-function(num) {
  if(num==1) {
    return(1)
  }
  return(2^(num-1))
}
neglogtrans<-function(num,lowest) {
  tymin = -1*log(lowest,2)
  return(-1*(-tymin+-1*log(num,2))/tymin)
}
xscale<-function(num,fac) {
  if(num==0) { 
    return(0)
  }
  return(num**(fac))
}

fac = 4
maxval = max(dg[,4])
#par(las=2)
par(mar=c(5,5,1,1))
#par(mfrow=c(2,1))
plot(1,type="n",ylim=c(1,logtrans(maxval*2)),xlim=c(0,1),yaxt='n',xaxt='n',xlab="fraction of expressed feature",ylab="Depth",cex.lab=1.5)
yvals = seq(0,logtrans(maxval*2),1)
axis(2,at=yvals,labels=lapply(yvals,untrans))
xvals = seq(0,1,0.01)
axis(1,at=lapply(xvals,xscale,fac),labels=xvals)
d = dg[order(dg[,4]),]
tot = sum(d[,3]-d[,2])
lines((cumsum(d[,3]-d[,2])/tot)**fac,lapply(d[,4],logtrans),col="#777777",lwd=4)
d = dinter[order(dinter[,4]),]
tot = sum(d[,3]-d[,2])
lines((cumsum(d[,3]-d[,2])/tot)**fac,lapply(d[,4],logtrans),col="#0000FFAA",lwd=4)
d = dintron[order(dintron[,4]),]
tot = sum(d[,3]-d[,2])
lines((cumsum(d[,3]-d[,2])/tot)**fac,lapply(d[,4],logtrans),col="#D04FD0AA",lwd=4)
d = dexon[order(dexon[,4]),]
tot = sum(d[,3]-d[,2])
lines((cumsum(d[,3]-d[,2])/tot)**fac,lapply(d[,4],logtrans),col="#FF0000AA",lwd=4)
#fac = 1
#plot(1,type="n",ylim=c(1,logtrans(maxval*2)),xlim=c(0,1),yaxt='n',xaxt='n',xlab="fraction of coveraged genome",ylab="Depth",cex.lab=1.5)
#yvals = seq(0,logtrans(maxval*2),1)
#axis(2,at=yvals,labels=lapply(yvals,untrans))
#xvals = seq(0,1,0.01)
#axis(1,at=lapply(xvals,xscale,fac),labels=xvals)
#d = dg[order(dg[,4]),]
#tot = sum(d[,3]-d[,2])
#lines((cumsum(d[,3]-d[,2])/tot)**fac,lapply(d[,4],logtrans),col="#777777",lwd=4)
#d = dinter[order(dinter[,4]),]
#tot2 = sum(d[,3]-d[,2])
#dif = tot-tot2
#lines(((cumsum(d[,3]-d[,2])+dif)/tot)**fac,lapply(d[,4],logtrans),col="#0000FFAA",lwd=4)
#d = dintron[order(dintron[,4]),]
#tot2 = sum(d[,3]-d[,2])
#dif = tot-tot2
#lines(((cumsum(d[,3]-d[,2])+dif)/tot)**fac,lapply(d[,4],logtrans),col="#D04FD0AA",lwd=4)
#d = dexon[order(dexon[,4]),]
#tot2 = sum(d[,3]-d[,2])
#dif = tot-tot2
#lines(((cumsum(d[,3]-d[,2])+dif)/tot)**fac,lapply(d[,4],logtrans),col="#FF0000AA",lwd=4)
dev.off()
