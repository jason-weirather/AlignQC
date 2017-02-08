#!/usr/bin/env Rscript

##### BEGIN FUNCTION BLOCK #####
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

  recwid = input_width
  axcex = text_adjust
  longest = 4000

  dexon = d[d[,2]=="exon",4]
  dintron = d[d[,2]=="intron",4]
  dintergenic = d[d[,2]=="intergenic",4]

  par(oma=c(0.5,0.5,0.5,0.5))
  par(mar=c(4.5,5,2,1))
  mat=rbind(c(1,2,3),c(4,4,4))
  layout(mat,c(1.5,5,1.1),c(6,6))


  exon = length(dexon)
  intron = length(dintron)
  intergenic = length(dintergenic)
  tot = exon+intron+intergenic
  plot(1,type="n",xlim=c(-50,500),ylim=c(0,tot*1.1),xaxt='n',ylab="Number of Reads",bty="n",xlab="",yaxt="n",cex.lab=axcex,yaxs="i")
  axis(2,lwd=recwid,cex.axis=axcex)
  rect(0,0,500,exon,col="#777777",lwd=recwid)
  rect(0,exon,500,exon+intron,col="#FFDB58",lwd=recwid)
  rect(0,exon+intron,500,exon+intron+intergenic,col="#007FFF",lwd=recwid)
  mtext("All",side=1,at=-100,adj=0,line=1,cex=axcex*0.67)

  # Do by length
  #find biggest
  biggest = 0
  for(i in seq(0,longest-500,500)) {
    currsize = length(d[d[,4]>i & d[,4]<=i+500,1])
    if(currsize > biggest) { biggest = currsize }
  }
  currsize = length(d[d[,4]>longest,1])
  if(currsize > biggest) { biggest = currsize }

  par(mar=c(4.5,2,2,1))
  plot(1,type="n",xlim=c(0,longest),ylim=c(0,biggest*1.1),ylab="",bty="n",xlab="Aligned Length (bp)",xaxt='n',yaxt='n',cex.lab=axcex,yaxs="i")
  axis(side=1,at=seq(0,longest,500),cex.axis=axcex,lwd=recwid)
  axis(side=2,lwd=recwid,cex.axis=axcex)
  for(i in seq(0,longest-500,500)) {
    exon = length(d[d[,2]=="exon" & d[,4]>i & d[,4]<=i+500,4])
    intron = length(d[d[,2]=="intron" & d[,4]>i & d[,4]<=i+500,4])
    intergenic = length(d[d[,2]=="intergenic" & d[,4]>i & d[,4]<=i+500,4])
    rect(i,0,i+500,exon,col="#777777",lwd=recwid)
    rect(i,exon,i+500,exon+intron,col="#FFDB58",lwd=recwid)
    rect(i,exon+intron,i+500,exon+intron+intergenic,col="#007FFF",lwd=recwid)
  }

  exon = length(d[d[,2]=="exon" & d[,4]>longest,1])
  intron = length(d[d[,2]=="intron" & d[,4]>longest,1])
  intergenic = length(d[d[,2]=="intergenic" & d[,4]>longest,1])
  tot = exon+intron+intergenic
  if(tot>0) {
    plot(1,type="n",xlim=c(-50,500),ylim=c(0,tot*1.1),ylab="",bty="n",xlab="",xaxt='n',yaxt="n",cex.lab=axcex,yaxs="i")
    axis(2,lwd=recwid,cex.axis=axcex)
    rect(0,0,500,exon,col="#777777",lwd=recwid)
    rect(0,exon,500,exon+intron,col="#FFDB58",lwd=recwid)
    rect(0,exon+intron,500,exon+intron+intergenic,col="#007FFF",lwd=recwid)
    mtext(paste(">",longest),side=1,at=500,adj=1,line=1,cex=axcex*0.67)
  } else {
    plot.new()
  }

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
    currsize = length(d[d[,3]==i &d[,2]=="exon",1])
    if(currsize > biggest) { biggest = currsize }
    currsize = length(d[d[,3]==i &d[,2]=="intron",1])
    if(currsize > biggest) { biggest = currsize }
    currsize = length(d[d[,3]==i &d[,2]=="intergenic",1])
    if(currsize > biggest) { biggest = currsize }
  }
  currsize = length(d[d[,3]>=maxexon & d[,2]=="exon",1])
  if(currsize > biggest) { biggest = currsize }
  currsize = length(d[d[,3]>=maxexon & d[,2]=="intron",1])
  if(currsize > biggest) { biggest = currsize }
  currsize = length(d[d[,3]>=maxexon & d[,2]=="intergenic",1])
  if(currsize > biggest) { biggest = currsize }

  endspace = 4
  par(mar=c(5,5,1.5,1))
  plot(1,type="n",xlim=c(0,maxexon+endspace),ylim=c(0,logtrans(biggest*2)),ylab="Number of Reads",bty="n",xlab="Exon Count",xaxt='n',cex.lab=axcex,yaxt='n',yaxs="i")
  axispoints = seq(0,logtrans(biggest*2),1)
  axis(2,at=axispoints,labels=lapply(axispoints,untrans),cex.axis=axcex,lwd=recwid)
  axis(1,at=seq(1,maxexon,1),labels=seq(1,maxexon,1),cex.axis=axcex,lwd=recwid)
  mtext(paste(">",maxexon),side=1,at=maxexon+endspace-1,line=1,cex=axcex*0.67)
  for(i in seq(1,maxexon)) {
    exon = length(d[d[,2]=="exon" & d[,3]==i,4])
    intron = length(d[d[,2]=="intron" & d[,3]==i,4])
    intergenic = length(d[d[,2]=="intergenic" & d[,3]==i,4])
    rect(i+0.1-0.5,0,i+0.8-0.5,logtrans(exon),col="#777777",lwd=recwid)
    rect(i+0-0.5,0,i+0.5-0.5,logtrans(intron),col="#FFDB58",lwd=recwid)
    rect(i+0.3-0.5,0,i+0.6-0.5,logtrans(intergenic),col="#007FFF",lwd=recwid)
  }
  # get the longest pooled
  exon = length(d[d[,2]=="exon" & d[,3]>maxexon,4])
  intron = length(d[d[,2]=="intron" & d[,3]>maxexon,4])
  intergenic = length(d[d[,2]=="intergenic" & d[,3]>maxexon,4])
  rect(maxexon+endspace-1+0.1-0.5,0,maxexon+endspace-1+0.8-0.5,logtrans(exon),col="#777777",lwd=recwid)
  rect(maxexon+endspace-1+0-0.5,0,maxexon+endspace-1+0.5-0.5,logtrans(intron),col="#FFDB58",lwd=recwid)
  rect(maxexon+endspace-1+0.3-0.5,0,maxexon+endspace-1+0.6-0.5,logtrans(intergenic),col="#007FFF",lwd=recwid)
  dev.off()
}
####### END FUNCTION BLOCK ##########


args=commandArgs(trailingOnly=TRUE)
if(length(args)<2) {
  stop("Must supply input and output\n",call.=FALSE)
}
file_name = args[2]
infile = args[1]
infilex = substr(infile,nchar(infile)-1,nchar(infile))
if(infilex=="gz") {
  d<-read.table(infile)
} else {
  d<-read.table(gzfile(infile))
}
input_width = 3
if(length(args)>2) { 
  input_width = args[3]
}
text_adjust = 1.75
if(length(args)>3) {
  text_adjust = args[4]
}
make_image(d,file_name,input_width,text_adjust)
