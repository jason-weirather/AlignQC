#!/usr/bin/env Rscript

### BEGIN FUNCTION BLOCK
make_image<-function(d,outfile,input_width,text_adjust) {
  # decide output type
  filex = substr(outfile,nchar(outfile)-2,nchar(outfile))
  if(filex=="pdf") {
    pdf(outfile,bg="#FFFFFF",width=12,height=7)
  } else if (filex=="png") {
    png(outfile,bg="#FFFFFF",width=823,height=480,units="px")
  } else {
      stop("Unsupported type for output file.\n",call.=FALSE)
  }

  par(mar=c(8,6,6,1))
  layout(matrix(c(1,2),,ncol=2),widths=c(1,4),heights=c(1,1))


  #plot correct vs not correct
  y0 = sum(d[d[,1]==0,2])
  y1 = sum(d[d[,1]!=0,2])
  yh= sum(abs(d[d[,1]>0,2]))

  yl = sum(abs(d[d[,1]<0,2]))
  print(d[d[,1]<0,2])
  print(yl)
  print(d[d[,1]>0,2])
  print(yh)
  ymax = max(y0,y1)
  plot(1,type='n',xlim=c(0,2),ylim=c(0,ymax),xaxt='n',xlab="",ylab="Number of exon boundries",bty="n",cex.lab=text_adjust,yaxt="n")
  if(y0 >0) { rect(0.25,0,1,y0,lwd=input_width,col="#86C67C") }
  if(yl > 0) {rect(1,0,1.75,yl,lwd=input_width,col="#B0E2FF") }
  if(yh >0) {rect(1,yl,1.75,yh+yl,lwd=input_width,col="#8968CD") }
  axis(side=1,at=c(0.5,1.5),labels=c('perfect','imperfect'),lwd=input_width,cex.axis=text_adjust,las=2)
  axis(side=2,lwd=input_width,cex.axis=text_adjust)

  par(mar=c(8,2,6,3))

  wrongvals = d[d[,1]!=0,]
  ymax = max(wrongvals[,2])
  xmin = min(wrongvals[,1])
  xmax = max(wrongvals[,1])
  plot(1,type='n',xlim=c(xmin,xmax),ylim=c(0,ymax),xaxt='n',xlab="distance (bp)",ylab="",bty="n",cex.lab=text_adjust,yaxt='n',)
  for(i in xmin:xmax) {
    step = 0.5
    mycol = "#B0E2FF"
    if(i!=0) {
      if(i>0) { 
        step=0.5; 
        mycol = "#8968CD"
      }
      if(length(wrongvals[wrongvals[,1]==i,2]) > 0 && wrongvals[wrongvals[,1]==i,2]!=0) {
        rect(i+step,0,i+1+step,wrongvals[wrongvals[,1]==i,2],lwd=input_width, col=mycol)
      }
    }
  }
  axlabs = c()
  if(xmin < 0) {
    for (i in seq(xmin,-1,1)) {
      axlabs = c(axlabs,i)
    }
  }
  if(xmax > 0) {
    for (i in seq(1,xmax,1)) {
      axlabs = c(axlabs,i)
    }
  }
  axis(side=1,at=seq(xmin,xmax,1),labels=seq(xmin,xmax,1),lwd=input_width,cex.axis=text_adjust)
  axis(side=2,lwd=input_width,cex.axis=text_adjust)
  mtext("truncated Exons",at=xmin,adj=0,side=1,line=3,cex=text_adjust)
  mtext("elongated Exons",at=xmax,adj=1,side=1,line=3,cex=text_adjust)
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

input_width = 3
if(length(args) > 2) {
  input_width = args[3]
}
text_adjust = 1.75
if(length(args) > 3) {
  text_adjust = args[4]
}
make_image(d1,outfile,input_width,text_adjust)
