#!/opt/R/3.2.1/bin/Rscript

# Read the inputs
# 1. bed depth file
# 2. reference lengths <chr> <length> per line tsv
# 3. output file (.png or .pdf)
args=commandArgs(trailingOnly=TRUE)
if(length(args)<3) {
  stop("Must supply input depth, reference lengths, and output files\n",call.=FALSE)
}

# decide output type
filex = substr(args[3],nchar(args[3])-2,nchar(args[3]))
if(filex=="pdf") {
  pdf(args[3],bg="#FFFFFF")
} else if (filex=="png") {
  png(args[3],bg="#FFFFFF")
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
g = unique(sort(e[,4]))
colfunc<-colorRampPalette(c("#00000088","#0000FF88","#FF000088"))
#colfunc<-colorRampPalette(c("#000000","#000033","#000066","#000099","#0000CC","#0000FF11","#FFCC00AA","#FF0000"))
colarray<-colfunc(length(f[,1]))
alldep = sort(rep(f[,3]-f[,2],f[,4]))
colong = colfunc(length(rep(f[,3]-f[,2],f[,4])))
gcol = vector(length=length(g))
for(i in 1:length(g)) {
  dind = match(g[i],alldep)
  gcol[i] = colong[dind]
  print(gcol)
}
cind = i
for(i in 1:length(f[,1])) {
#for(i in 1:25000) {
  cind = match(f[i,1],ordered_names)
  clen = ordered_lens[cind]
  fcol = gcol[f[i,4]]
  rect(cind-0.4,mytrans(f[i,2],clen,minlen,maxlen),cind+0.4,mytrans(f[i,3],clen,minlen,maxlen),col=fcol,border=fcol,lwd=0.01)
}
par(mar=c(10,0,10,5))
plot(1,type="n",xlim=c(0,2),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
mindepth = min(g)
maxdepth = max(g)
z = 0
step = 1/length(f[,4])
rect(1-0.4,0,1+0.4,1)
mult = 100
for(i in seq(1,length(f[,4]),mult)) {
  rect(1-0.4,z,1+0.4,z+step,col=gcol[i],border=gcol[f[i,4]])
  z = step*mult + z
}
ylabs = vector(length=length(g))
yval = vector(length=length(g))
for(i in 1:length(f[,4])) {
  gi = match(f[i,4],g,1)
  ylabs[gi] = f[i,4]
  yval[gi] = i/length(f[,4])
}
prev = 0
for(i in 1:length(ylabs)) {
  if(yval[i]>prev+0.1 & yval[i]<0.95) {
    mtext(ylabs[i],side=4,at=yval[i])
    prev = yval[i]
  }
}
mtext(max(ylabs),side=4,at=1)
curr = 1
#axis(4,at=yval,labels=ylabs)
dev.off()
