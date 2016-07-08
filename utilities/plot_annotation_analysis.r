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
  png(args[2],bg="#FFFFFF")
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
longest = 4000

full = length(d[d[,1]=="full",1])
partial = length(d[d[,1]=="partial",1])
unannot = length(d[d[,1]=="unannotated",1])

par(mar=c(4,4,1.5,0.5),oma=c(0.5,1,0.5,1))
mat=rbind(c(1,2,3),c(4,4,4))
layout(mat,c(1,5,1),c(1,1))


tot = full+partial+unannot
plot(1,type="n",xlim=c(0,500),ylim=c(0,tot*1.1),xaxt='n',ylab="Counts any length",bty="n",xlab="",cex.axis=1.2,cex.lab=1.2)
rect(0,0,500,partial,col="#777777",lwd=recwid)
rect(0,partial,500,full+partial,col="#FF0000",lwd=recwid)
rect(0,full+partial,500,full+partial+unannot,col="#C6E2FF",lwd=recwid)
mtext("All reads",side=1,at=-100,adj=0)

# Do by length
#find biggest
biggest = 0
for(i in seq(0,longest-500,500)) {
  currsize = length(d[d[,2]>i & d[,2]<=i+500,1])
  if(currsize > biggest) { biggest = currsize }
}
currsize = length(d[d[,2]>longest,1])
if(currsize > biggest) { biggest = currsize }

plot(1,type="n",xlim=c(0,longest),ylim=c(0,biggest*1.1),ylab="Count reads per aligned length",bty="n",xlab="Aligned length (bp)",xaxt='n',cex.axis=1.2,cex.lab=1.2)
axis(side=1,at=seq(0,longest,500),cex.axis=1.2)
for(i in seq(0,longest-500,500)) {
  full = length(d[d[,1]=="full" & d[,2]>i & d[,2]<=i+500,1])
  partial = length(d[d[,1]=="partial" & d[,2]>i & d[,2]<=i+500,1])
  unannot = length(d[d[,1]=="unannotated" & d[,2]>i & d[,2]<=i+500,1])
  rect(i,0,i+500,partial,col="#777777",lwd=recwid)
  rect(i,partial,i+500,full+partial,col="#FF0000",lwd=recwid)
  rect(i,full+partial,i+500,full+partial+unannot,col="#C6E2FF",lwd=recwid)
}

full = length(d[d[,1]=="full" & d[,2]>longest,1])
partial = length(d[d[,1]=="partial" & d[,2]>longest,1])
unannot = length(d[d[,1]=="unannotated" & d[,2]>longest,1])
tot = full+partial+unannot
plot(1,type="n",xlim=c(0,500),ylim=c(0,tot*1.1),ylab="Count longest reads",bty="n",xlab="",xaxt='n',cex.axis=1.2,cex.lab=1.2)
rect(0,0,500,partial,col="#777777",lwd=recwid)
rect(0,partial,500,full+partial,col="#FF0000",lwd=recwid)
rect(0,full+partial,500,full+partial+unannot,col="#C6E2FF",lwd=recwid)
mtext(paste(">",longest),side=1,at=-100,adj=0)


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
maxexon = 30
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

recwid = 1
endspace = 4
plot(1,type="n",xlim=c(0,maxexon+endspace),ylim=c(0,logtrans(biggest*2)),ylab="Count reads",bty="n",xlab="Exon count (alignment)",xaxt='n',cex.axis=1.2,cex.lab=1.2,yaxt='n')
axispoints = seq(0,logtrans(biggest*2),1)
axis(2,at=axispoints,labels=lapply(axispoints,untrans),cex.axis=1.2)
axis(1,at=seq(1,maxexon,1),labels=seq(1,maxexon,1),cex.axis=1.2)
mtext(paste(">",maxexon),side=1,at=maxexon+endspace-1)
for(i in seq(1,maxexon)) {
  full = length(d[d[,1]=="full" & d[,3]==i,1])
  partial = length(d[d[,1]=="partial" & d[,3]==i,1])
  unannot = length(d[d[,1]=="unannotated" & d[,3]==i,1])
  rect(i+0.1-0.5,0,i+0.8-0.5,logtrans(partial),col="#777777",lwd=recwid)
  rect(i+0-0.5,0,i+0.5-0.5,logtrans(full),col="#FF0000",lwd=recwid)
  rect(i+0.3-0.5,0,i+0.6-0.5,logtrans(unannot),col="#C6E2FF",lwd=recwid)
}
# get the longest pooled
full = length(d[d[,1]=="full" & d[,3]>maxexon,1])
partial = length(d[d[,1]=="partial" & d[,3]>maxexon,1])
unannot = length(d[d[,1]=="unannotated" & d[,3]>maxexon,1])
rect(maxexon+endspace-1+0.1-0.5,0,maxexon+endspace-1+0.8-0.5,logtrans(partial),col="#777777",lwd=recwid)
rect(maxexon+endspace-1+0-0.5,0,maxexon+endspace-1+0.5-0.5,logtrans(full),col="#FF0000",lwd=recwid)
rect(maxexon+endspace-1+0.3-0.5,0,maxexon+endspace-1+0.6-0.5,logtrans(unannot),col="#C6E2FF",lwd=recwid)
