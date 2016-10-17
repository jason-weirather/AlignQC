#!/usr/bin/env Rscript

##### BEGIN FUNCTION BLOCK ######
make_image<-function(d,outfile,input_width,text_adjust) {
  # decide output type
  filex = substr(outfile,nchar(outfile)-2,nchar(outfile))
  if(filex=="pdf") {
    pdf(outfile,bg="#FFFFFF")
  } else if (filex=="png") {
    png(outfile,bg="#FFFFFF")
  } else {
      stop("Unsupported type for output file.\n",call.=FALSE)
  }

  pcol = '#777777'
  fcol = '#7FC97F'
  ucol = '#FDC086'

  recwid = input_width
  axcex = text_adjust
  longest = 4000

  full = length(d[d[,1]=="full",1])
  partial = length(d[d[,1]=="partial",1])
  unannot = length(d[d[,1]=="unannotated",1])

  par(oma=c(0.5,0.5,0.5,0.5))
  par(mar=c(4,5,2,1))
  mat=rbind(c(1,2,3),c(4,4,4))
  layout(mat,c(1.5,5,1.1),c(6,6))


  tot = full+partial+unannot
  plot(1,type="n",xlim=c(-50,500),ylim=c(0,tot*1.1),xaxt='n',ylab="Number of Reads",bty="n",xlab="",yaxt='n',yaxs='i',cex.lab=axcex)
  axis(2,lwd=recwid,cex.axis=axcex)
  rect(0,0,500,partial,col=pcol,lwd=recwid)
  rect(0,partial,500,full+partial,col=fcol,lwd=recwid)
  rect(0,full+partial,500,full+partial+unannot,col=ucol,lwd=recwid)
  mtext("All",side=1,at=-100,adj=0,line=1,cex=axcex*0.67)

  # Do by length
  #find biggest
  biggest = 0
  for(i in seq(0,longest-500,500)) {
    currsize = length(d[d[,2]>i & d[,2]<=i+500,1])
    if(currsize > biggest) { biggest = currsize }
  }
  currsize = length(d[d[,2]>longest,1])
  if(currsize > biggest) { biggest = currsize }

  par(mar=c(4,2,2,1))
  plot(1,type="n",xlim=c(0,longest),ylim=c(0,biggest*1.1),ylab="",bty="n",xlab="Reference Transcript Length (bp)",xaxt='n',yaxt='n',yaxs='i',cex.lab=axcex)
  axis(side=1,at=seq(0,longest,500),cex.axis=axcex,lwd=recwid)
  axis(side=2,lwd=recwid,cex.axis=axcex)
  for(i in seq(0,longest-500,500)) {
    full = length(d[d[,1]=="full" & d[,2]>i & d[,2]<=i+500,1])
    partial = length(d[d[,1]=="partial" & d[,2]>i & d[,2]<=i+500,1])
    unannot = length(d[d[,1]=="unannotated" & d[,2]>i & d[,2]<=i+500,1])
    rect(i,0,i+500,partial,col=pcol,lwd=recwid)
    rect(i,partial,i+500,full+partial,col=fcol,lwd=recwid)
    rect(i,full+partial,i+500,full+partial+unannot,col=ucol,lwd=recwid)
  }

  full = length(d[d[,1]=="full" & d[,2]>longest,1])
  partial = length(d[d[,1]=="partial" & d[,2]>longest,1])
  unannot = length(d[d[,1]=="unannotated" & d[,2]>longest,1])
  tot = full+partial+unannot
  plot(1,type="n",xlim=c(-50,500),ylim=c(0,tot*1.1),ylab="Count longest reads",bty="n",xlab="",xaxt='n',yaxt='n',yaxs='i',cex.lab=axcex)
  axis(2,lwd=recwid,cex.axis=axcex)
  rect(0,0,500,partial,col=pcol,lwd=recwid)
  rect(0,partial,500,full+partial,col=fcol,lwd=recwid)
  rect(0,full+partial,500,full+partial+unannot,col=ucol,lwd=recwid)
  mtext(paste(">",longest),side=1,at=500,adj=1,line=1,cex=axcex*0.67)


  logtrans<-function(num) {
    if(num==0) {
      return(0)
    }
    if(num==1) {
      return(1)
    }
    return(log(num,2)+1)
  }
  untrans<-function(num) {
    if(num==0) {
      return(0)
    }
    if(num==1) {
      return(1)
    }
    return(2^(num-1))
  }
  neglogtrans<-function(num,lowest) {
    tymin = -1*log(lowest,2)
    return(-1*(-tymin+-1*log(num,2))/tymin)
  }


  # plot by exon distribution
  #find biggest
  maxexon = 20
  biggest = 0
  for(i in seq(1,maxexon,1)) {
    currsize = length(d[d[,3]==i &d[,1]=="full",1])
    if(currsize > biggest) { biggest = currsize }
    currsize = length(d[d[,3]==i &d[,1]=="partial",1])
    if(currsize > biggest) { biggest = currsize }
    currsize = length(d[d[,3]==i &d[,1]=="unannotated",1])
    if(currsize > biggest) { biggest = currsize }
  }
  currsize = length(d[d[,3]>=maxexon & d[,1]=="full",1])
  if(currsize > biggest) { biggest = currsize }
  currsize = length(d[d[,3]>=maxexon & d[,1]=="partial",1])
  if(currsize > biggest) { biggest = currsize }
  currsize = length(d[d[,3]>=maxexon & d[,1]=="unannotated",1])
  if(currsize > biggest) { biggest = currsize }

  endspace = 4
  par(mar=c(5,5,1.5,1))
  plot(1,type="n",xlim=c(0,maxexon+endspace),ylim=c(0,logtrans(biggest*2)),ylab="Number of Reads",bty="n",xlab="Exon Count (Reference Transcript)",xaxt='n',cex.lab=axcex,yaxt='n',yaxs='i')
  axispoints = seq(0,logtrans(biggest*2),1)
  axis(2,at=axispoints,labels=lapply(axispoints,untrans),cex.axis=axcex,lwd=recwid)
  axis(1,at=seq(1,maxexon,1),labels=seq(1,maxexon,1),cex.axis=axcex,lwd=recwid)
  mtext(paste(">",maxexon),side=1,at=maxexon+endspace-1,line=1,cex=axcex*0.67)
  for(i in seq(1,maxexon)) {
    full = length(d[d[,1]=="full" & d[,3]==i,1])
    partial = length(d[d[,1]=="partial" & d[,3]==i,1])
    unannot = length(d[d[,1]=="unannotated" & d[,3]==i,1])
    rect(i+0.1-0.5,0,i+0.8-0.5,logtrans(partial),col=pcol,lwd=recwid)
    rect(i+0-0.5,0,i+0.5-0.5,logtrans(full),col=fcol,lwd=recwid)
    rect(i+0.3-0.5,0,i+0.6-0.5,logtrans(unannot),col=ucol,lwd=recwid)
  }
  # get the longest pooled
  full = length(d[d[,1]=="full" & d[,3]>maxexon,1])
  partial = length(d[d[,1]=="partial" & d[,3]>maxexon,1])
  unannot = length(d[d[,1]=="unannotated" & d[,3]>maxexon,1])
  rect(maxexon+endspace-1+0.1-0.5,0,maxexon+endspace-1+0.8-0.5,logtrans(partial),col=pcol,lwd=recwid)
  rect(maxexon+endspace-1+0-0.5,0,maxexon+endspace-1+0.5-0.5,logtrans(full),col=fcol,lwd=recwid)
  rect(maxexon+endspace-1+0.3-0.5,0,maxexon+endspace-1+0.6-0.5,logtrans(unannot),col=ucol,lwd=recwid)
  dev.off()
}

args=commandArgs(trailingOnly=TRUE)
if(length(args)<2) {
  stop("Must supply input and output\n",call.=FALSE)
}
outfile=args[2]
infile = args[1]
infilex = substr(infile,nchar(infile)-1,nchar(infile))
if(infilex=="gz") {
  d1<-read.table(infile)
} else {
  d1<-read.table(gzfile(infile))
}
d<-data.frame(as.character(d1[,5]))
d<-cbind(d,as.numeric(d1[,12]),as.numeric(d1[,9]))
input_width = 3
if(length(args) > 2) {
  input_width = args[3]
}
text_adjust = 1.75
if(length(args) > 3) {
  text_adjust = args[4]
}
make_image(d,outfile,input_width,text_adjust)
