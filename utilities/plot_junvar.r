#!/opt/R/3.2.1/bin/Rscript

### BEGIN FUNCTION BLOCK
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

  par(mar=c(8,6,6,1))
  layout(matrix(c(1,2),,ncol=2),widths=c(1.5,3),heights=c(1,1))


  #plot correct vs not correct
  y0 = sum(d[d[,1]==0,2])
  y1 = sum(d[d[,1]!=0,2])
  ymax = max(y0,y1)
  plot(1,type='n',xlim=c(0,2),ylim=c(0,ymax),xaxt='n',xlab="",ylab="Number of exon boundries",bty="n",cex.lab=text_adjust,yaxt="n")
  rect(0.25,0,1,y0,lwd=input_width,col="#86C67C")
  rect(1,0,1.75,y1,lwd=input_width,col="#B0E2FF")
  axis(side=1,at=c(0.5,1.5),labels=c('perfect','imperfect'),lwd=input_width,cex.axis=text_adjust,las=2)
  axis(side=2,lwd=input_width,cex.axis=text_adjust)

  par(mar=c(8,2,6,2))

  wrongvals = d[d[,1]!=0,]
  ymax = max(wrongvals[,2])
  xmax = max(wrongvals[,1])
  plot(1,type='n',xlim=c(1,xmax),ylim=c(0,ymax),xaxt='n',xlab="distance (bp)",ylab="",bty="n",cex.lab=text_adjust,yaxt='n',)
  for(i in 1:xmax) {
    rect(i,0,i+1,wrongvals[i,2],lwd=input_width, col="#B0E2FF")
  }
  axis(side=1,at=wrongvals[,1],labels=wrongvals[,1],lwd=input_width,cex.axis=text_adjust)
  axis(side=2,lwd=input_width,cex.axis=text_adjust)
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
