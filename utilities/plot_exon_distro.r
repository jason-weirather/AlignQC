#!/usr/bin/env Rscript

######## BEGIN FUNCTION BLOCK ##########
make_image<-function(d,outfile,input_width,text_adjust) {
  # decide output type
  filex = substr(outfile,nchar(outfile)-2,nchar(outfile))
  if(filex=="pdf") {
    pdf(outfile,bg="#FFFFFF",height=4.5,width=7)
  } else if (filex=="png") {
    png(outfile,bg="#FFFFFF",width=480,height=240)
  } else {
    stop("Unsupported type for output file.\n",call.=FALSE)
  }


  recwid = input_width
  axcex = text_adjust

  par(oma=c(0.5,0.5,0.5,0.5))
  par(mar=c(4,5,2,0))


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
    currsize = length(d[d[,2]==i,1])
    if(currsize > biggest) { biggest = currsize }
  }
  currsize = length(d[d[,2]>=maxexon,1])
  if(currsize > biggest) { biggest = currsize }

  endspace = 4
  plot(1,type="n",xlim=c(0,maxexon+endspace),ylim=c(0,logtrans(biggest*2)+1),ylab="Number of Reads",bty="n",xlab="Exon Count",xaxt='n',yaxt='n',yaxs='i',cex.lab=axcex,yaxt='n')
  axispoints = seq(0,logtrans(biggest*2),1)
  axis(2,at=axispoints,labels=lapply(axispoints,untrans),cex.axis=axcex,lwd=recwid)
  axis(1,at=seq(1,maxexon,1),labels=seq(1,maxexon,1),cex.axis=axcex,lwd=recwid)
  mtext(paste(">",maxexon),side=1,at=maxexon+endspace-1.5,line=1,cex=axcex)
  for(i in seq(1,maxexon)) {
    exon = length(d[d[,2]==i,1])
    rect(i+0.1-0.5,0,i+0.8-0.5,logtrans(exon),col="#777777",lwd=recwid)
  }
  # get the longest pooled
  exon = length(d[d[,2]>maxexon,1])
  rect(maxexon+endspace-1+0.1-0.5,0,maxexon+endspace-1+0.8-0.5,logtrans(exon),col="#777777",lwd=recwid)
  dev.off()
}

args=commandArgs(trailingOnly=TRUE)
if(length(args)<2) {
  stop("Must supply input and output\n",call.=FALSE)
}
infile = args[1]
outfile = args[2]
infilex = substr(infile,nchar(infile)-1,nchar(infile))
if(infilex=="gz") {
  d<-read.table(infile)
} else {
  d<-read.table(gzfile(infile))
}
input_width = 3
if(length(args) > 2) {
  input_width = args[3]
}
text_adjust = 1.25
if(length(args) > 3) {
  text_adjust = args[4]
}
make_image(d,outfile,input_width,text_adjust)
