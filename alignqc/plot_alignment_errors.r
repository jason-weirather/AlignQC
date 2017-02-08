#!/usr/bin/env Rscript


args=commandArgs(trailingOnly=TRUE)
#args=c('report13.txt','')
if(length(args)<2) {
  stop("Must supply input and output files\n",call.=FALSE)
}

filex = substr(args[2],nchar(args[2])-2,nchar(args[2]))

rins = c()
rmis = c()
rdel = c()
if (length(args) > 2) {
  if(length(args) != 8) {
    stop("If ranges are defined, must define 6 after input and output file.\n",call.=FALSE)
  }
  rins = c(as.numeric(args[3]),as.numeric(args[4]))
  rmis = c(as.numeric(args[5]),as.numeric(args[6]))
  rdel = c(as.numeric(args[7]),as.numeric(args[8]))
}

if(filex=="pdf") {
  pdf(args[2])
} else if (filex=="png") {
  png(args[2])
} else {
    stop("Unsupported type for output file.\n",call.=FALSE)
}
d<-read.table(args[1],header=TRUE)
mat = matrix(seq(1,30,1),nrow=5,ncol=6,byrow=TRUE)
#par(mfrow=c(5,6))
layout(mat,widths=c(1,1,1,1,1,1),heights=c(1,1,1,1,1,1))
par(bg="#FFFFFF")
par(pty="m")
par(mar=c(0.1,0.1,0.1,0.1))
par(oma=c(2,8,8,4))
basecex = 3
scalecex = 2.5
baselabcex = 1.8
scalelabcex=1.8

baselabcol="#555555"
scalelabcol="#555555"

mismatches = d[which(d$target != '-' & d$query != '-' &
               d$target != d$query),][,3]/d$total[1]
if(length(rmis)==0) {
  rmis = c(min(mismatches),max(mismatches))
}
if(min(mismatches)==max(mismatches)){
  rmis = c(max(0,min(mismatches)-0.000001),min(1,max(mismatches)+0.00001))
}
rmispal = colorRampPalette(c("blue","white","red"))(100)

ins = d[which(d$target == '-' &
               d$target != d$query),][,3]/d$total[1]
if(length(rins)==0) {
  rins = c(min(ins),max(ins))
}
if(min(ins)==max(ins)){
  rins = c(max(0,min(ins)-0.000001),min(1,max(ins)+0.00001))
}

rinspal = colorRampPalette(c("#7570B3","#FFFFFF","#E7298A"))(100)
del = d[which(d$query == '-' &
               d$target != d$query),][,3]/d$total[1]
if(length(rdel)==0) {
  rdel = c(min(del),max(del))
}
if(min(del)==max(del)){
  rdel = c(max(0,min(del)-0.000001),min(1,max(del)+0.00001))
}
rdelpal = colorRampPalette(c("#1B9E77","#FFFFFF","#D95F02"))(100)

a = c('-','A','C','G','T')
for (i in 1:length(a)) {
  for (j in 1:length(a)) {
    plot(1,type="n",xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',bty='n',ylab="",xlab="")
    if(i==5 & j==1) {
      mtext('Reference',side=2,line=5,cex=baselabcex,adj=0,at=0,col=baselabcol)
    }
    if(i==1) {
      if(a[j]!='-') {
        mtext(a[j],side=3,line=1,cex=basecex)
      } else {
        mtext("del",side=3,line=1,cex=basecex)
      }
    }
    if(j==1) {
      if(a[i]!='-') {
        mtext(a[i],side=2,line=1,cex=basecex,las=1)
      } else {
        mtext("ins",side=2,line=1,cex=basecex,las=1)
      }
    }
    if(i==1 & j==5) {
      mtext('Query',side=3,line=5,cex=baselabcex,at=1,adj=1,col=baselabcol)
    }
    num = d[which(d$target==a[i] & d$query==a[j]),3]/d$total[1]
    vcol = "#000000"
    if(a[i]=='-') {
      vcol = rinspal[max(1,100*(num-rins[1])/(rins[2]-rins[1]))]
    } else if (a[j] == '-') {
      vcol = rdelpal[max(1,100*(num-rdel[1])/(rdel[2]-rdel[1]))]
    } else {
      vcol = rmispal[max(1,100*(num-rmis[1])/(rmis[2]-rmis[1]))]
    }
    if (i!=j) {
      rect(0,0,1,1,col=vcol)
    }
  }
  if(i==1) {
    plot(1,type="n",xlim=c(0,1),ylim=c(rins[1],rins[2]),xaxt='n',yaxt='n',bty='n',ylab="",xlab="")
    anums = seq(rins[1],rins[2],(rins[2]-rins[1])/2)
    for(i in 1:length(anums)) {
      anums[i] = round(anums[i],digits=4)
    }
    axis(side=4,pos=0.2,cex.axis=scalecex,las=1,at=anums,lwd.ticks=5)
    legend_image =as.raster(matrix(rev(rinspal)),ncol=1)
    rasterImage(legend_image,0,rins[1],0.2,rins[2])
    mtext('insertion',side=3,adj=0,line=0.8,cex=scalelabcex,col=scalelabcol)  
  } else if(i==3) {
    plot(1,type="n",xlim=c(0,1),ylim=c(rmis[1],rmis[2]),xaxt='n',yaxt='n',bty='n',ylab="",xlab="")
    anums = seq(rmis[1],rmis[2],(rmis[2]-rmis[1])/2)
    for(i in 1:length(anums)) {
      anums[i] = round(anums[i],digits=4)
    }
    axis(side=4,pos=0.2,cex.axis=scalecex,las=1,at=anums,lwd.ticks=5)
    legend_image =as.raster(matrix(rev(rmispal)),ncol=1)
    rasterImage(legend_image,0,rmis[1],0.2,rmis[2])    
    mtext('mismatch',side=3,adj=0,line=0.8,cex=scalelabcex,col=scalelabcol)  
  } else if(i==5) {
    plot(1,type="n",xlim=c(0,1),ylim=c(rdel[1],rdel[2]),xaxt='n',yaxt='n',bty='n',ylab="",xlab="")
    anums = seq(rdel[1],rdel[2],(rdel[2]-rdel[1])/2)
    for(i in 1:length(anums)) {
      anums[i] = round(anums[i],digits=4)
    }
    axis(side=4,pos=0.2,cex.axis=scalecex,las=1,at=anums,lwd.ticks=5)
    legend_image =as.raster(matrix(rev(rdelpal)),ncol=1)
    rasterImage(legend_image,0,rdel[1],0.2,rdel[2])  
    mtext('deletion',side=3,adj=0,line=0.8,cex=scalelabcex,col=scalelabcol)  
  } else {
    plot(1,type="n",xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',bty='n',ylab="",xlab="")
  }
}
dev.off()
