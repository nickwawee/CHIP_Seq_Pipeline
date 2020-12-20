#Common Regions
#This script will find common regions between heatmaps and output a bedfile that contains all common regions. It also outputs a new matrix for each binned matrix that contains a common bed region
#Its input are the bedfiles from each sample as well as the matrix file for each sample
rm(list=ls())
#libraries
library(parallel)
#functions
bednames<-function(bedmat){#this function takes a BED3 file and returns a single character vector of names for regions
  bedmat[,1]=paste(bedmat[,1],rep(":", nrow(bedmat)),sep="")
  bedmat[,2]=paste(as.character(bedmat[,2]), as.character(bedmat[,3]),sep="-")
  bednames<-paste(bedmat[,1],bedmat[,2],sep="")      
}
hm="H3K27ac"
dt="seqdepthnorm"
#paths
machine="server"
if (machine=="local"){
  p="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/bedfiles/"#bedfiles path
  p2="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/matrixfiles/"
  nc=4
}
if (machine=="server"){
  p=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",hm,"_",dt, "/bedfiles/", sep="")#bedfiles path
  p2=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",hm,"_", dt,"/matrixfiles/", sep="")
  nc=28
}
#parameters
refpoint="TSS"
#loading
bedfiles<-list.files(path=p, pattern=paste("sortedregions_refpoint_",refpoint,sep=""))
bedlist<-mclapply(1:length(bedfiles), function(i){
  bednames(read.table(file=paste(p, bedfiles[i], sep=""))[,1:3])
}, mc.cores=nc)

listlengths<-rep(NA, length(bedlist))
for (i in 1:length(bedlist)){
  listlengths[i]<-length(bedlist[[i]])
}
if (length(duplicated(listlengths))==length(listlengths)){
  commonbed=bedlist[[1]]
}else{commonbed<-bedlist[[which(listlengths==min(listlengths))]]}
if (length(which(duplicated(commonbed)))>0){
  commonbed<-commonbed[-which(duplicated(commonbed))]
}

#loading matrix files
matfiles=list.files(path = p2, pattern=refpoint)
matfiles=matfiles[grep(pattern=".tab", x=matfiles)]
if (length(grep(pattern="common", x=matfiles))>0){
  matfiles=matfiles[-grep(pattern="common", x=matfiles)]
}
mclapply(1:length(bedfiles), function(i){
  oldmatrix=read.csv(paste(p2,matfiles[i], sep=""), sep="\t",header=FALSE, stringsAsFactors = FALSE, skip=3)
  newmatrix<-oldmatrix[match(commonbed, bedlist[[i]]),]
  row.names(newmatrix)<-commonbed
  write.table(newmatrix, paste(p2,"common_",matfiles[i],sep=""), row.names = T, sep="\t")
}, mc.cores = nc, mc.silent = T)

write.table(commonbed, paste(p,"",refpoint,"_commonbed.bed", sep=""), quote=F)

#TSS
#parameters
refpoint="TSS"
#loading
bedfiles<-list.files(path=p, pattern=paste("sortedregions_refpoint_",refpoint,sep=""))
bedlist<-mclapply(1:length(bedfiles), function(i){
  bednames(read.table(file=paste(p, bedfiles[i], sep=""))[,1:3])
}, mc.cores=nc)

listlengths<-rep(NA, length(bedlist))
for (i in 1:length(bedlist)){
  listlengths[i]<-length(bedlist[[i]])
}
if (length(duplicated(listlengths))==length(listlengths)){
  commonbed=bedlist[[1]]
}else{commonbed<-bedlist[[which(listlengths==min(listlengths))]]}
if (length(which(duplicated(commonbed)))>0){
  commonbed<-commonbed[-which(duplicated(commonbed))]
}

#loading matrix files
matfiles=list.files(path = p2, pattern=refpoint)
matfiles=matfiles[grep(pattern=".tab", x=matfiles)]
if (length(grep(pattern="common", x=matfiles))>0){
  matfiles=matfiles[-grep(pattern="common", x=matfiles)]
}
mclapply(1:length(bedfiles), function(i){
  oldmatrix=read.csv(paste(p2,matfiles[i], sep=""), sep="\t",header=FALSE, stringsAsFactors = FALSE, skip=3)
  newmatrix<-oldmatrix[match(commonbed, bedlist[[i]]),]
  row.names(newmatrix)<-commonbed
  write.table(newmatrix, paste(p2,"common_",matfiles[i],sep=""), row.names = T, sep="\t")
}, mc.cores = nc, mc.silent = T)

write.table(commonbed, paste(p,"",refpoint,"_commonbed.bed", sep=""), quote=F)
~                                                                              
