#!/usr/bin/env Rscript

# Read the inputs
# 1. line_plot_table.txt.gz
# 2. total_distro_table.txt.gz
# 3. chr_distro_table.txt.gz
# 4. output (pdf or png)
args=commandArgs(trailingOnly=TRUE)
if(length(args)<4) {
  stop("Must supply input \n",call.=FALSE)
}
# decide output type
filex = substr(args[4],nchar(args[4])-2,nchar(args[4]))
if(filex=="pdf") {
  pdf(args[4],bg="#FFFFFF")
} else if (filex=="png") {
  png(args[4],bg="#FFFFFF")
} else {
    stop("Unsupported type for output file.\n",call.=FALSE)
}



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

layout(rbind(c(1,2),c(3,4),c(5,5)),widths=c(1.25,4),heights=c(1,1,1))
par(las=1)

# get lowest y value for any chrom
#d<-read.table(gzfile(args[3]))
absolute_min = 0.00001
#ymin = min(d[,4]/d[,5])

d<-read.table(gzfile(args[2]))
par(mar=c(4,9,1,0.5))
####### first plot the coverage
ymax = max(d[,3]/d[,4])
#ymin = min(d[,3]/d[,4])
ylocs = c(0.000001,0.00001,0.0001,0.001,0.01,0.1,1)
ytlocs = lapply(ylocs,neglogtrans,absolute_min)
plot(1,type="n",xlim=c(1,1),ylim=c(0,1),xaxt="n",xlab="",ylab="",yaxt="n",bty="n")
axis(2,at=ytlocs,labels=ylocs,col="#FF0000",col.axis="#BB0000",line=3)
axis(2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
par(las=0)
mtext("Aligned Fraction",side=2,line=7)
par(las=1)

cov = d[,3]
total = d[,4]
height = cov[1]/total[1]
theight = neglogtrans(max(height,absolute_min),absolute_min)
rect(1-0.4,0,1,theight,col="#FF0000")
rect(1,0,1+0.4,height,col="#777777")
par(las=2)
mtext("all",1,at=1,line=0.5,adj=1)
par(las=1)


par(mar=c(4,0.5,1,1))
####### first plot the coverage
d<-read.table(gzfile(args[3]))
chr_list = unique(d[,1])
ymax = max(d[,4]/d[,5])
#absolute_min = 0.0001
ymin = min(d[,4]/d[,5])
ylocs = c(0.000001,0.00001,0.0001,0.001,0.01,0.1,1)
ytlocs = lapply(ylocs,neglogtrans,absolute_min)
plot(1,type="n",xlim=c(1,length(chr_list)),ylim=c(0,1),xaxt="n",xlab="",ylab="",yaxt="n",bty="n")
#axis(2,at=ytlocs,labels=ylocs,col="#FF0000",col.axis="#BB0000",line=3)
#axis(2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
par(las=0)
#mtext("Aligned Fraction",side=2,line=6)
par(las=1)
z = 0
for (chr in chr_list) {
  z = z+1
  cov = d[which(d[,1]==chr),4]
  total = d[which(d[,1]==chr),5]
  height = cov[1]/total[1]
  theight = neglogtrans(max(height,absolute_min),absolute_min)
  rect(z-0.4,0,z,theight,col="#FF0000")
  rect(z,0,z+0.4,height,col="#777777")
  par(las=2)
  mtext(chr,1,at=z,line=0.5,adj=1)
  par(las=1)
}

##### plot box depth for all data
d<-read.table(gzfile(args[2]))
par(mar=c(4,8,1,0.5))
ymax = max(d[,1])
ymaxtrans = logtrans(ymax)
plot(1,type="n",xlim=c(1,1),ylim=c(1,ymaxtrans),xaxt="n",yaxt="n",xlab="",ylab="Depth",bty="n",cex.lab=1.5)
ylocs = seq(1,ymaxtrans,1)
ylabs=lapply(ylocs,untrans)
axis(2,labels=ylabs,at=ylocs)
depths = lapply(d[,1],logtrans)
cnts = d[,2]
boxplot(unlist(rep(depths,cnts)),at=1,add=TRUE,outline=FALSE,yaxt="n",frame=FALSE)
par(las=2)
mtext("all",1,at=1,line=0.5,adj=1)
par(las=1)


d<-read.table(gzfile(args[3]))
ymax = max(d[,2])
ymaxtrans = logtrans(ymax)
par(mar=c(4,0.5,1,1))
plot(1,type="n",xlim=c(1,length(chr_list)),ylim=c(1,ymaxtrans),xaxt="n",yaxt="n",xlab="",ylab="Depth",bty="n")
ylocs = seq(1,ymaxtrans,1)
ylabs = lapply(ylocs,untrans)
#axis(2,labels=ylabs,at=ylocs)
z = 0
for (chr in chr_list) {
  z = z+1
  #depths = lapply(lapply(d[which(d[,1]==chr),2],logtrans),round)
  #depths = d[which(d[,1]==chr),2]
  depths = lapply(d[which(d[,1]==chr),2],logtrans)
  cnts = d[which(d[,1]==chr),3]
  boxplot(unlist(rep(depths,cnts)),at=z,add=TRUE,outline=FALSE,yaxt="n",frame=FALSE)
  par(las=2)
  mtext(chr,1,at=z,line=0.5,adj=1)
  par(las=1)
}

d<-read.table(gzfile(args[1]))
par(mar=c(5,8,1,0.5))
ymax=max(d[,1])
ymaxtrans=logtrans(ymax)
xmin=min(d[,2]/d[,3])
plot(1,type="n",ylim=c(1,ymaxtrans),xlim=c(0,max(d[,2]/d[,3])-xmin),yaxt='n',xlab="fraction of genome",ylab="Depth",cex.lab=1.5)
lines(d[,2]/d[,3]-xmin,lapply(d[,1],logtrans),lwd=3)
ylocs = seq(1,ymaxtrans,1)
ylabs = lapply(ylocs,untrans)
axis(2,at=ylocs,labels=ylabs)
dev.off()
