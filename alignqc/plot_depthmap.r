#!/usr/bin/env Rscript

# Read the inputs
# 1. bed depth file
# 2. reference lengths <chr> <length> per line tsv
# 3. strata key
# 4. output file (.png or .pdf)
args=commandArgs(trailingOnly=TRUE)
if(length(args)<4) {
  stop("Must supply input depth, reference lengths, strata key, and output files\n",call.=FALSE)
}

# decide output type
filex = substr(args[4],nchar(args[4])-2,nchar(args[4]))
if(filex=="pdf") {
  pdf(args[4],bg="#FFFFFF")
} else if (filex=="png") {
  png(args[4],bg="#FFFFFF")
} else {
    stop("Unsupported type for output file.\n",call.=FALSE)
}



mytrans<-function(pos,clen,minlen,maxlen) {
  # find range it should cover
  tot = max((9/10)*(clen/maxlen)+(1/10),(1/10))
  return((pos/clen)*tot)
}

d<-read.table(args[2])
#ordered_names = as.vector(d[order(names(chrlens)),1])
#ordered_lens = as.vector(d[order(names(chrlens)),2])
ordered_names = as.vector(d[order(d[,1]),1])
ordered_lens = as.vector(d[order(d[,1]),2])
minlen = min(ordered_lens)
maxlen = max(ordered_lens)
chr = cbind.data.frame(ordered_names,ordered_lens)
longest = max(chr[,2])

dstrata<-read.table(args[3],header=TRUE)

layout(rbind(c(1,2)),heights=c(1),widths=c(4,1))
par(mar=c(5,1,1,1))
plot(1,type="n",xlim=c(1,length(ordered_names)),ylim=c(0,1),bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
par(las=2)
mtext(ordered_names,side=1,at=1:length(ordered_names))
par(las=1)
for (i in 1:length(ordered_names)) {
  clen = ordered_lens[i]
  rect(i-0.4,0,i+0.4,mytrans(chr[i,2],clen,minlen,maxlen))
}
e<-read.table(gzfile(args[1]))
f<-e[order(e[,4]),]
colfunc<-colorRampPalette(c("#00000088","#0000FF88","#FF000088"))
colarray<-colfunc(length(dstrata[,1]))
cind = i
for(i in 1:length(f[,1])) {
  cind = match(f[i,1],ordered_names)
  clen = ordered_lens[cind]
  #print(f[i,4])
  fcol = colarray[f[i,4]]
  #print(fcol)
  rect(cind-0.4,mytrans(f[i,2],clen,minlen,maxlen),cind+0.4,mytrans(f[i,3],clen,minlen,maxlen),col=fcol,border=fcol,lwd=0.01)
}
par(mar=c(10,0,10,5))
plot(1,type="n",xlim=c(0,2),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
mindepth = min(dstrata[1,])
maxdepth = max(dstrata[1,])
z = 0
step = 1/length(f[,4])
legend(-1,1,legend=dstrata[,1],fill=rev(colarray),xpd=TRUE,bty='n')
dev.off()
