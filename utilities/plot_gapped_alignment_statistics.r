#!/usr/bin/env Rscript

#### BEGIN FUNCTION DEFINITIONS ####
make_image<-function(d,infile,title_string,input_width,text_adjust) {
  filex = substr(infile,nchar(infile)-2,nchar(infile))
  if(filex=="pdf") {
    pdf(infile,bg="#FFFFFF")
  } else if(filex=="png") {
    png(infile,bg="#FFFFFF")
  } else {
    stop("Unsupported type of output file\n",call.=FALSE)
  }  

  #par(bg="#FFFFFF99")
  par(oma=c(0.5,0.5,0.5,0.5))
  par(mar=c(2,5,2,1))
  mat=rbind(c(1,2,3),c(4,5,6))
  layout(mat,c(1.5,5,1.1),c(6,6))

  recwid = input_width
  axcex = text_adjust
  single = length(d[d[,2]=="original",1])
  gapped = length(d[d[,2]=="gapped",1])
  transchimera = length(d[d[,2]=="chimera",1])
  selfchimera = length(d[d[,2]=="self-chimera" | d[,2]=="self-chimera-atypical",1])
  unaligned = length(d[d[,2]=="unaligned",1])
  tot = gapped+transchimera+selfchimera+single+unaligned
  acols = c("#FFFFFF","#777777","#FF0000")
  plot(1,type="n",xlim=c(-50,500),ylim=c(0,tot*1.1),xaxt='n',ylab="Number of Reads",bty="n",xlab="",yaxt="n",cex.lab=axcex,yaxs="i")
  axis(2,lwd=recwid,cex.axis=axcex)
  rect(0,0,500,single,col="#777777",lwd=recwid)
  rect(0,single,500,single+gapped,col="#FF0000",lwd=recwid)
  rect(0,single+gapped,500,single+gapped+selfchimera,col="#551A8B",lwd=recwid)
  rect(0,single+gapped+selfchimera,500,single+gapped+selfchimera+transchimera,col="#007FFF",lwd=recwid)
  rect(0,single+gapped+selfchimera+transchimera,500,single+gapped+selfchimera+transchimera+unaligned,col="#FFFFFF",lwd=recwid)
  mtext("All",side=1,at=-100,adj=0,line=1,cex=axcex*0.67)

  #mtext(title_string,side=3,outer=TRUE,line=-2,cex=axcex)
  ###### Read length-wise bar plot for split alignmetns #####
  longest = 4000

  #find the biggest bin
  biggest = 0
  for(i in seq(0,longest-500,500)) {
    currsize = length(d[d[,5]>i & d[,5]<=i+500,1])
    if(currsize > biggest) { biggest = currsize }
  }
  currsize = length(d[d[,5]>longest,1])
  if(currsize > biggest) { biggest = currsize }


  par(mar=c(2,2,2,1))
  plot(1,type="n",xlim=c(0,longest),ylim=c(0,biggest*1.1),ylab="",bty="n",xlab="",xaxt='n',yaxt='n',cex.lab=axcex,yaxs="i")
  axis(side=1,at=seq(0,longest,500),cex.axis=axcex,lwd=recwid)
  axis(2,lwd=recwid,cex.axis=axcex)
  for(i in seq(0,longest-500,500)) {
    single = length(d[d[,2]=="original" & d[,5]>i & d[,5]<=i+500,1])
    rect(i,0,i+500,single,col="#777777",lwd=recwid)
    gapped = length(d[d[,2]=="gapped" & d[,5]>i & d[,5]<=i+500,1])
    rect(i,single,i+500,single+gapped,col="#FF0000",lwd=recwid)
    selfchimera = length(d[(d[,2]=="self-chimera" | d[,2]=="self-chimera-atypical")  & d[,5]>i & d[,5]<=i+500,1])
    rect(i,single+gapped,i+500,single+gapped+selfchimera,col="#551A8B",lwd=recwid)
    transchimera = length(d[d[,2]=="chimera" & d[,5]>i & d[,5]<=i+500,1])
    rect(i,single+gapped+selfchimera,i+500,single+gapped+selfchimera+transchimera,col="#007FFF",lwd=recwid)
    unalign = length(d[d[,2]=="unaligned" & d[,5]>i & d[,5]<=i+500,1])
    rect(i,single+gapped+selfchimera+transchimera,i+500,single+gapped+selfchimera+transchimera+unalign,col="#FFFFFF",lwd=recwid)
  }

  ###last case###
  single = length(d[d[,2]=="original" & d[,5]>longest,1])
  gapped = length(d[d[,2]=="gapped" & d[,5]>longest,1])
  selfchimera = length(d[(d[,2]=="self-chimera" | d[,2]=="self-chimera-atypical") & d[,5]>longest,1])
  transchimera = length(d[d[,2]=="chimera" & d[,5]>longest,1])
  unalign = length(d[d[,2]=="unaligned" & d[,5]>longest,1])
  tot = single+gapped+selfchimera+transchimera+unalign
  if(tot>0) {
    plot(1,type="n",xlim=c(-50,500),ylim=c(0,tot*1.1),ylab="",bty="n",xlab="",xaxt='n',yaxt="n",cex.lab=axcex,yaxs="i")
    axis(2,lwd=recwid,cex.axis=axcex)
    rect(0,0,500,single,col="#777777",lwd=recwid)
    rect(0,single,500,single+gapped,col="#FF0000",lwd=recwid)
    rect(0,single+gapped,500,single+gapped+selfchimera,col="#551A8B",lwd=recwid)
    rect(0,single+gapped+selfchimera,500,single+gapped+selfchimera+transchimera,col="#007FFF",lwd=recwid)
    rect(0,single+gapped+selfchimera+transchimera,500,single+gapped+selfchimera+transchimera+unalign,col="#FFFFFF",lwd=recwid)
    mtext(paste(">",longest),side=1,at=500,adj=1,line=1,cex=axcex*0.67)
  } else {
    plot.new()
  }

  par(mar=c(5,5,1.5,1))
  plot(1,type='n',xlim=c(-50,500),ylim=c(-0.04,1.04),bty="n",xaxt='n',ylab="Fraction of Read Aligned",xlab="",yaxt="n",yaxs='i',cex.lab=axcex)
  axis(2,lwd=recwid,cex.axis=axcex)
  dat1 = d[d[,2]!="unaligned",4]
  dat2 = d[d[,2]!="unaligned",5]
  dcom = dat1/dat2
  boxplot(dcom,add=TRUE,at=250,boxwex=900,border="#000000",frame=FALSE,axes=FALSE,outline=FALSE,boxlwd=recwid,whisklty=1,whisklwd=recwid,staplelwd=recwid)
  mtext("All",side=1,at=0,adj=0,line=1,cex=axcex*0.67)


  ####### Plot the distribution of reads ###########
  par(mar=c(5,2,1.5,1))
  plot(1,type='n',xlim=c(0,longest),ylim=c(-0.04,1.04),bty="n",xaxt='n',ylab="",xlab="Read length(bp)",yaxt='n',cex.lab=axcex,yaxs="i")
  axis(side=1,at=seq(0,longest,500),cex.axis=axcex,lwd=recwid)
  axis(2,lwd=recwid,cex.axis=axcex)
  for(i in seq(0,longest-500,500)) {
    dat1 = d[d[,2]!="unaligned" & d[,5]>i & d[,5]<=i+500,4]
    dat2 = d[d[,2]!="unaligned" & d[,5]>i & d[,5]<=i+500,5]
    dcom = dat1/dat2
    boxplot(dcom,add=TRUE,at=i+250,boxwex=900,border="#000000",frame=FALSE,axes=FALSE,outline=FALSE,whisklty=1,boxlwd=recwid,staplelwd=recwid,whisklwd=recwid)
  }
  dat1 = d[d[,2]!="unaligned" & d[,5]>longest,4]
  if(length(dat1)>0){
    dat2 = d[d[,2]!="unaligned" & d[,5]>longest,5]
    dcom = dat1/dat2
    plot(1,type='n',xlim=c(-50,500),ylim=c(-0.04,1.04),bty="n",xaxt='n',ylab="",xlab="",yaxt="n",cex.lab=axcex,yaxs='i')
    axis(2,lwd=recwid,cex.axis=axcex)
    boxplot(dcom,add=TRUE,at=250,boxwex=900,border="#000000",frame=FALSE,axes=FALSE,outline=FALSE,col="#00000000",boxlwd=recwid,whisklty=1,whisklwd=recwid,staplelwd=recwid)
    mtext(paste(">",longest),side=1,at=500,adj=1,line=1,cex=axcex*0.67)
  } else {
    plot.new()
  }
  dev.off()
}
####### END FUNCTION DEFINITIONS #########

args=commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Must supply input and output files\n <input_file (.gz supported)> <output_file (.pdf or .png)> <line width (default 3)> <text R cex adjust (default 1.75)>\n",call.=FALSE)
}
file_name = args[2]
infilex = substr(args[1],nchar(args[1])-1,nchar(args[1]))
if(infilex=="gz") {
  d<-read.table(args[1])
} else {
  d<-read.table(gzfile(args[1]))
}

line_width = 3
text_adjust = 1.75
if(length(args)>2) { 
  line_width=args[3]
}
if(length(args)>3) {
  text_adjust = args[4]
}

make_image(d,file_name,"",line_width,text_adjust)
