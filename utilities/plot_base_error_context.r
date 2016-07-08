#!/opt/R/3.2.1/bin/Rscript

args=commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Must supply input and output files\n",call.=FALSE)
}

filex = substr(args[2],nchar(args[2])-2,nchar(args[2]))

rins = c()
rmis = c()
rdel = c()
if (length(args) > 2) {
  if(length(args) != 8) {
    stop("If ranges are defined, must define 6 after input and output file.  rinsmin rinsmax rmismin rmismax rdelmin rdelmax\n",call.=FALSE)
  }
  rins = c(as.numeric(args[3]),as.numeric(args[4]))
  rmis = c(as.numeric(args[5]),as.numeric(args[6]))
  rdel = c(as.numeric(args[7]),as.numeric(args[8]))
}  
if(filex=="pdf") {
  pdf(args[2])
} else if(filex=="png") {
  png(args[2],width=8,height=8,res=300,units="in")
} else {
    stop("unsupported type of output file.  rinsmin rinsmax rmismin rmismax rdelmin rdelmax\n",call.=FALSE)
}
d<-read.table(args[1],header=TRUE)
par(bg="#FFFFFF")
par(mfrow=c(5,6))
par(mar=c(0.1,0.1,0.1,0.1))
par(oma=c(2,8,8,4))


contextcex= 1.5
contextlabcex = 1.2
basecex = 2
baselabcex = 1.8
scalecex=2.5
scalelabcex=1.8

baselabcol="#555555"
contextlabcol="#555555"
scalelabcol="#555555"

### Find the ranges for base mismatches ###
mismatches = d[which(d$reference != '-' & d$query != '-' &
               d$reference != d$query),][,5]
if(length(rmis)==0) {
  rmis = c(min(mismatches),max(mismatches))
}
if(min(mismatches)==max(mismatches)){
  rmis = c(max(0,min(mismatches)-0.000001),min(1,max(mismatches)+0.00001))
}

rmispal = colorRampPalette(c("blue","white","red"))(100)

### Find the ranges for insertions ###
ins = d[which(d$reference == '-' &
               d$reference != d$query),][,5]
if(length(rins)==0) {
  rins = c(min(ins),max(ins))
}
if(min(ins)==max(ins)){
  rins = c(max(0,min(ins)-0.000001),min(1,max(ins)+0.00001))
}

rinspal = colorRampPalette(c("#7570B3","#FFFFFF","#E7298A"))(100)

### Find the ranges for deletions ###
del = d[which(d$query == '-' &
               d$reference != d$query),][,5]
if(length(rdel)==0) {
  rdel = c(min(del),max(del))
}
if(min(del)==max(del)){
  rdel = c(max(0,min(del)-0.000001),min(1,max(del)+0.00001))
}
rdelpal = colorRampPalette(c("#1B9E77","#FFFFFF","#D95F02"))(100)

types = c('-','A','C','G','T')
afterbases = c('A','C','G','T')
beforebases = c('T','G','C','A')
for (ci in 1:5) { # the reference 
  for (cj in 1:5) { # what the query became
    if (types[ci] != types[cj]) {
      plot(1,type="n",xlim=c(0,4),ylim=c(0,4),xaxt='n',yaxt='n',ann=FALSE,bty="n")
      for (i in 1:4) { # the before base
        for (j in 1:4) { # the after base
          v = d[which(d$reference==types[ci] & d$query==types[cj] & 
                      d$before==beforebases[i] & d$after==afterbases[j]),]
          num = v[1,5]
          vcol = "#000000"
          if(types[ci]=='-') {
            vcol = rinspal[max(1,100*(num-rins[1])/(rins[2]-rins[1]))]
          } else if (types[cj]=='-') {          
            vcol = rdelpal[max(1,100*(num-rdel[1])/(rdel[2]-rdel[1]))]
          } else { 
            vcol = rmispal[max(1,100*(num-rmis[1])/(rmis[2]-rmis[1]))]
          }
          rect(j-1,i-1,j,i,col=vcol)
        }
      }
      if (ci==1 | cj+1 ==ci) {
        mtext(c('A','C','G','T'),side=3,at=seq(0.5,3.5,1),cex=contextcex,line=-0.4)
      }
      if (cj==1 | ci+1 ==cj) {
        mtext(rev(c('A','C','G','T')),side=2,at=seq(0.5,3.5,1),las=1,cex=contextcex,line=-0.1)
      }
    } else {
      plot(1,type="n",xlim=c(0,4),ylim=c(0,4),xaxt='n',yaxt='n',ann=FALSE,bty="n")
      if(ci==1 & cj==1) {
        mtext("before",side=2,adj=0,cex=contextlabcex,col=contextlabcol)
        mtext("after",side=3,adj=1,cex=contextlabcex,col=contextlabcol)
      }
    }
    if(ci==1) {
      if(types[cj]!='-') {
        mtext(types[cj],side=3,line=2,cex=basecex)
      } else {
        mtext("del",side=3,line=2,cex=basecex)
      }
    }
    if(cj==1) {
      if(types[ci]!='-') {
        mtext(types[ci],side=2,line=2,las=1,cex=basecex)
      } else {
        mtext("ins",side=2,line=2,las=1,cex=basecex)
      }
    }
    if(ci==5 & cj==1) {
      mtext('Reference',side=2,line=5.0,cex=baselabcex,at=0,adj=0,col=baselabcol)
    }
    if(ci==1 & cj==5) {
      mtext('Query',side=3,line=5.0,cex=baselabcex,at=4,adj=1,col=baselabcol)
    }
  }
  # the last entry
  ### Add legends
  if(ci==1) {
    plot(1,type='n',xlim=c(0,1),ylim=c(rins[1],rins[2]),xaxt='n',yaxt='n',ann=FALSE,bty="n",xaxs="i")
    anums = seq(rins[1],rins[2],(rins[2]-rins[1])/2)
    for(i in 1:length(anums)) {
      anums[i] = round(anums[i],digits=4)
    }
    axis(side=4,pos=0.2,cex.axis=scalecex,las=1,at=anums,lwd.ticks=5)
    legend_image = as.raster(matrix(rev(rinspal)),ncol=1)
    rasterImage(legend_image,0,rins[1],0.2,rins[2])
    mtext('insertion',side=3,adj=0,line=0.8,cex=scalelabcex,col=scalelabcol)
  }
  else if(ci==3) {
    plot(1,type='n',xlim=c(0,1),ylim=c(rmis[1],rmis[2]),xaxt='n',yaxt='n',ann=FALSE,bty="n")
    anums = seq(rmis[1],rmis[2],(rmis[2]-rmis[1])/2)
    for(i in 1:length(anums)) {
      anums[i] = round(anums[i],digits=4)
    }
    axis(side=4,pos=0.2,cex.axis=scalecex,las=1,at=anums,lwd.ticks=5)
    legend_image = as.raster(matrix(rev(rmispal)),ncol=1)
    rasterImage(legend_image,0,rmis[1],0.2,rmis[2])
    mtext('mismatch',side=3,adj=0,line=0.8,cex=scalelabcex,col=scalelabcol)
  }
  else if(ci==5) {
    plot(1,type='n',xlim=c(0,1),ylim=c(rdel[1],rdel[2]),xaxt='n',yaxt='n',ann=FALSE,bty="n")
    anums = seq(rdel[1],rdel[2],(rdel[2]-rdel[1])/2)
    for(i in 1:length(anums)) {
      anums[i] = round(anums[i],digits=4)
    }
    axis(side=4,pos=0.2,cex.axis=scalecex,las=1,at=anums,lwd.ticks=5)
    legend_image = as.raster(matrix(rev(rdelpal)),ncol=1)
    rasterImage(legend_image,0,rdel[1],0.2,rdel[2])
    mtext('deletion',side=3,adj=0,line=0.8,cex=scalelabcex,col=scalelabcol)
  }
  else {
    plot(1,type='n',xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',ann=FALSE,bty="n")
  }
}
dev.off()
