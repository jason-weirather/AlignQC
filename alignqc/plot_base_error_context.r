#!/usr/bin/env Rscript

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
  png(args[2],height=7,width=7,units="in",res=300)
} else {
    stop("unsupported type of output file.  rinsmin rinsmax rmismin rmismax rdelmin rdelmax\n",call.=FALSE)
}
d<-read.table(args[1],header=TRUE)
par(bg="#FFFFFF")
a = 0.5 # width of box
f = 1.5 # width of legend
b = 3 # height of box
e = 0.5 # height of padding row
par(mar=c(0.1,0.1,0.1,0.1))
par(oma=c(2,8,8,1))
layout(rbind(
 c(1,2,3,4,25),     # Insertions followed by insertion legend
 c(28,28,28,28,28), # Padding row
 c(5,6,7,8,26),     # Deletions followed by deletion legend
 c(29,29,29,29,29), # Padding row
 c(9,10,11,12,27),  # Mismatch row followed by mismatch legend 
 c(13,14,15,16,30), # Mismatch row followed by padding
 c(17,18,19,20,30),
 c(21,22,23,24,30)
),widths=c(
  a,a,a,a,f
),heights=c(
  b,e,b,e,b,b,b,b
))

contextoffsettop = -0.3
contextcex= 1.4
contextlabcex = 1.2
basecex = 2
baselabcex = 1.8
scalecex=2.5
scalelabcex=1.8

baselabcol="#333333"
contextlabcol="#444444"
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

types = c('A','C','G','T')
afterbases = c('A','C','G','T')
beforebases = c('T','G','C','A')
#Do insertions first
for (cj in 1:4) {
  plot(1,type="n",xlim=c(0,4),ylim=c(0,4),xaxt='n',yaxt='n',ann=FALSE,bty="n")
  for ( i in 1:4) { # before base
    for (j in 1:4) { # after base
      v = d[which(d$reference=='-' & d$query==types[cj] &
                  d$before==beforebases[i] & d$after==afterbases[j]),]
      num = v[1,5]
      vcol = "#000000"
      vcol = rinspal[max(1,100*(num-rins[1])/(rins[2]-rins[1]))]
      rect(j-1,i-1,j,i,col=vcol)
    }
  }
  mtext(c('A','C','G','T'),side=3,at=seq(0.5,3.5,1),cex=contextcex,line=contextoffsettop,col=contextlabcol)
  if (cj==1) {
    mtext(rev(c('A','C','G','T')),side=2,at=seq(0.5,3.5,1),las=1,cex=contextcex,line=-0.1,col=contextlabcol)
    mtext("insertion",side=2,line=4,cex=basecex,at=3)
    mtext("before",side=2,line=1.2,cex=contextcex,at=1.75,col=contextlabcol)
    mtext("after",side=3,line=1.1,cex=contextcex,at=1.3,col=contextlabcol)
  }
  if (cj==4) {
    mtext("observed error",side=3,at=-2,cex=basecex,line=4.5)
  }
  mtext(types[cj],side=3,line=2.2,cex=basecex)
}
#Do deletions next
for (ci in 1:4) {
  plot(1,type="n",xlim=c(0,4),ylim=c(0,4),xaxt='n',yaxt='n',ann=FALSE,bty="n")
  for ( i in 1:4) { # before base
    for (j in 1:4) { # after base
      v = d[which(d$reference==types[ci] & d$query=='-' &
                  d$before==beforebases[i] & d$after==afterbases[j]),]
      num = v[1,5]
      vcol = "#000000"
      vcol = rdelpal[max(1,100*(num-rdel[1])/(rdel[2]-rdel[1]))]
      rect(j-1,i-1,j,i,col=vcol)
    }
  }
  mtext(c('A','C','G','T'),side=3,at=seq(0.5,3.5,1),cex=contextcex,line=contextoffsettop,col=contextlabcol)
  if (ci==1) {
    mtext(rev(c('A','C','G','T')),side=2,at=seq(0.5,3.5,1),las=1,cex=contextcex,line=-0.1,col=contextlabcol)
    mtext("deleltion",side=2,line=4,cex=basecex,at=1)
  }
}
#Now do the mismatches
for (ci in 1:4) {
  for (cj in 1:4) {
    if (types[ci] == types[cj]) {
      # leave a blank plot and skip case they are equal
      plot(1,type="n",xlim=c(0,4),ylim=c(0,4),xaxt='n',yaxt='n',ann=FALSE,bty="n")
      if(cj == 1) {
        mtext(types[ci],side=2,line=2,las=1,cex=basecex)
      }
      next
    }
    # They are different ... a mismatch
    plot(1,type="n",xlim=c(0,4),ylim=c(0,4),xaxt='n',yaxt='n',ann=FALSE,bty="n")
    for (i in 1:4) { # The before base
      for (j in 1:4) { # The after base
        v = d[which(d$reference==types[ci] & d$query==types[cj] & 
                    d$before==beforebases[i] & d$after==afterbases[j]),]
        num = v[1,5]
        vcol = "#000000"
        vcol = rmispal[max(1,100*(num-rmis[1])/(rmis[2]-rmis[1]))]
        rect(j-1,i-1,j,i,col=vcol)
      }
    }
    if(ci==1 && cj>1) {
      mtext(c('A','C','G','T'),side=3,at=seq(0.5,3.5,1),cex=contextcex,line=contextoffsettop,col=contextlabcol)
    }
    if (ci>1 && cj==1) {
      mtext(rev(c('A','C','G','T')),side=2,at=seq(0.5,3.5,1),las=1,cex=contextcex,line=-0.1,col=contextlabcol)
    }
    if(cj == 1) {
      mtext(types[ci],side=2,line=2,las=1,cex=basecex)
    }
    if(ci==4 && cj==1) {
      mtext('mismatch',side=2,line=4,cex=basecex,at=4)
    }
  }
}
legendwidth = 0.075
# Now plot legends
par(mar=c(1,0,1,0))
plot(1,type='n',xlim=c(0,1),ylim=c(rins[1],rins[2]),xaxt='n',yaxt='n',ann=FALSE,bty="n")
anums = seq(rins[1],rins[2],(rins[2]-rins[1])/2)
for(i in 1:length(anums)) {
  anums[i] = round(anums[i],digits=4)
}
axis(side=4,pos=legendwidth,cex.axis=scalecex,las=1,at=anums,lwd.ticks=5)
legend_image = as.raster(matrix(rev(rinspal)),ncol=1)
rasterImage(legend_image,0,rins[1],legendwidth,rins[2])
#mtext('insertion',side=3,adj=0,line=0.8,cex=scalelabcex,col=scalelabcol)

plot(1,type='n',xlim=c(0,1),ylim=c(rdel[1],rdel[2]),xaxt='n',yaxt='n',ann=FALSE,bty="n")
anums = seq(rdel[1],rdel[2],(rdel[2]-rdel[1])/2)
for(i in 1:length(anums)) {
  anums[i] = round(anums[i],digits=4)
}
axis(side=4,pos=legendwidth,cex.axis=scalecex,las=1,at=anums,lwd.ticks=5)
legend_image = as.raster(matrix(rev(rdelpal)),ncol=1)
rasterImage(legend_image,0,rdel[1],legendwidth,rdel[2])
#mtext('deletion',side=3,adj=0,line=0.8,cex=scalelabcex,col=scalelabcol)

plot(1,type='n',xlim=c(0,1),ylim=c(rmis[1],rmis[2]),xaxt='n',yaxt='n',ann=FALSE,bty="n")
anums = seq(rmis[1],rmis[2],(rmis[2]-rmis[1])/2)
for(i in 1:length(anums)) {
  anums[i] = round(anums[i],digits=4)
}
axis(side=4,pos=legendwidth,cex.axis=scalecex,las=1,at=anums,lwd.ticks=5)
legend_image = as.raster(matrix(rev(rmispal)),ncol=1)
rasterImage(legend_image,0,rmis[1],legendwidth,rmis[2])
#mtext('mismatch',side=3,adj=0,line=0.8,cex=scalelabcex,col=scalelabcol)


dev.off()
