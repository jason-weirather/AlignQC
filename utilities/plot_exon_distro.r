#!/opt/R/3.3.0/bin/R

args=commandArgs(trailingOnly=TRUE)
if(length(args)<2) {
  stop("Must supply input and output\n",call.=FALSE)
}

# decide output type
filex = substr(args[2],nchar(args[2])-2,nchar(args[2]))
if(filex=="pdf") {
  pdf(args[2],bg="#FFFFFF")
} else if (filex=="png") {
  png(args[2],bg="#FFFFFF",width=480,height=240)
} else {
    stop("Unsupported type for output file.\n",call.=FALSE)
}

infilex = substr(args[1],nchar(args[1])-1,nchar(args[1]))
if(infilex=="gz") {
  d<-read.table(args[1])
} else {
  d<-read.table(gzfile(args[1]))
}

recwid = 1

par(mar=c(4,4,1.5,0.5),oma=c(0.5,1,0.5,1))


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
maxexon = 60
biggest = 0
for(i in seq(1,maxexon,1)) {
  currsize = length(d[d[,2]==i,1])
  if(currsize > biggest) { biggest = currsize }
}
currsize = length(d[d[,2]>=maxexon,1])
if(currsize > biggest) { biggest = currsize }

recwid = 1
endspace = 4
plot(1,type="n",xlim=c(0,maxexon+endspace),ylim=c(0,logtrans(biggest*2)),ylab="Count reads",bty="n",xlab="Exon count (alignment)",xaxt='n',cex.axis=1.2,cex.lab=1.2,yaxt='n')
axispoints = seq(0,logtrans(biggest*2),1)
axis(2,at=axispoints,labels=lapply(axispoints,untrans),cex.axis=1.2)
axis(1,at=seq(1,maxexon,1),labels=seq(1,maxexon,1),cex.axis=1.2)
mtext(paste(">",maxexon),side=1,at=maxexon+endspace-1)
for(i in seq(1,maxexon)) {
  exon = length(d[d[,2]==i,1])
  rect(i+0.1-0.5,0,i+0.8-0.5,logtrans(exon),col="#777777",lwd=recwid)
}
# get the longest pooled
exon = length(d[d[,2]>maxexon,1])
rect(maxexon+endspace-1+0.1-0.5,0,maxexon+endspace-1+0.8-0.5,logtrans(exon),col="#777777",lwd=recwid)
dev.off()
