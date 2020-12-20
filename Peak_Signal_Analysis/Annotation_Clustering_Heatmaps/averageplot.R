##Average plot
#This script will create a matrix and a line plot that takes the average values of a specific phenotype and will plot them based on the observed differentially regulated clusters
rm(list=ls())
#functions
bednames<-function(bedmat){#this function takes a BED3 file and returns a single character vector of names for regions
  bedmat[,1]=paste(bedmat[,1],rep(":", nrow(bedmat)),sep="")
  bedmat[,2]=paste(as.character(bedmat[,2]), as.character(bedmat[,3]),sep="-")
  bednames<-paste(bedmat[,1],bedmat[,2],sep="")      
}

#libraries
library(ggplot2)
library(reshape2)
library(parallel)

#parameters
ncores=40
before=2000
after=1000
histonemark="H3K27ac"
plot_ext=".png"
rt1="promoter"#region type.
machine="server"
dt="seqdepthnorm"

if (rt1=="promoter"){
  rp="TSS"
  rt2="Promoter"
}
if (rt1=="enhancer"){
  rp="center"#reference point
  rt2="Enhancer"
}
#paths
if (machine=="local"){
  p="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/matrixfiles/"
  p2="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/bedfiles/"
  p3="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/plots/"
  p4="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/mousedata.tab"
}
if (machine=="server"){
  p=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",histonemark,"_",dt,"/matrixfiles/",sep="")
  p2=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",histonemark,"_",dt,"/bedfiles/",sep="")
  p3=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",histonemark,"_",dt,"/plots/",sep="")
  p4="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/mousedata.tab"
}
#regions and cluster numbers, feature data
features<-read.table(p4, header=T, sep="\t")
bedfiles<-list.files(path=p2, pattern=rt1)

if (rt1=="enhancer"){
  allbedfiles<-mclapply(1:length(bedfiles),function(i){
    test<-cbind(bednames(read.csv(paste(p2,bedfiles[i],sep=""), header = F, sep="\t")), rep(strsplit(strsplit(bedfiles[i], split="r")[[1]][3], split=".b")[[1]][1]))
  }, mc.cores = ncores)
}
if (rt1=="promoter"){
  allbedfiles<-mclapply(1:length(bedfiles),function(i){
    test<-cbind(bednames(read.csv(paste(p2,bedfiles[i],sep=""), header = F, sep="\t")), rep(strsplit(strsplit(bedfiles[i], split="r")[[1]][4], split=".b")[[1]][1]))
  }, mc.cores = ncores)
}

clusterdf<-as.data.frame(do.call(rbind, allbedfiles))
colnames(clusterdf)<-c("names", "cluster")
clusterdf<-clusterdf[which(clusterdf$names%in%as.character(read.table(paste(p2,rp,"_commonbed.bed",sep=""))$x)),]

ptvec=c("Wildtype", "Lean Nnat", "Huge Nnat")#, "Huge NnatD9")
clustervec<-levels(factor(clusterdf$cluster))

newlist<-mclapply(1:length(clustervec), function(l){
  cluster=clustervec[l]
  newclusterdf<-clusterdf[which(clusterdf$cluster==cluster),]
  gene_num=as.character(nrow(newclusterdf))
  
  #matrices and bedfile
  matfiles=list.files(path = p, pattern=rp)
  matfiles=matfiles[grep(pattern=".tab", x=matfiles)]
  matfiles=matfiles[grep(pattern="common", x=matfiles)]

  combined<-mclapply(1:length(matfiles), function(i){
    filen=matfiles[i]
    n=paste(strsplit(filen, split="_")[[1]][4])#makes name 
    tb=read.table(paste(p,filen,sep=""),sep="\t",stringsAsFactors = FALSE)#reads matrix file
    tb<-tb[which(rownames(tb)%in%clusterdf$names[which(clusterdf$cluster==cluster)]),]
  }, mc.cores = ncores/2)
  names<-unlist(mclapply(1:length(matfiles), function(i){
    filen=matfiles[i]
    n=paste(strsplit(filen, split="_")[[1]][4])#makes name 
  }, mc.cores = ncores/2))
  names(combined)<-names
  combinedlist<-mclapply(1:length(ptvec), function(j){
    pt=ptvec[j]
    typenums<-features$mousenumber[which(features$genopheno==pt)]
    type<-combined[which(names(combined)%in%as.character(typenums))]#returns a list of dataframes that contains cluster regions and binnned values
    typedf<-do.call(rbind,type)#combines list into a matrix
    typeavglist<-mclapply(1:length(newclusterdf$names), function(i){
      matchname<-as.character(newclusterdf$names[i])
      grep(matchname, rownames(typedf))
      filtdf<-colMeans(typedf[grep(matchname, rownames(typedf)),])
      }, mc.cores = ncores/4)
    names(typeavglist)<-as.character(newclusterdf$names)
    typeavgdf<-do.call(rbind,typeavglist)#dataframe that contains an average region enrichment value for all replicates per cluster
    typeavgvec<-colMeans(typeavgdf, na.rm = T)
    df<-rbind(rep(pt, ncol(typeavgdf)),typeavgvec, seq(-before/1000, after/1000, length.out=length(typeavgvec)))
    row.names(df)<-c("Type", "AvgValue", "BinNum")
    return(t(df))
  }, mc.cores=ncores/2)
  combineddf<-as.data.frame(do.call(rbind, combinedlist))
  combineddf$AvgValue<-as.numeric(as.character(combineddf$AvgValue))
  combineddf$BinNum<-as.numeric(as.character(combineddf$BinNum))
  combineddf$cluster<-rep(paste("Cluster ",cluster,", [",gene_num,"]",sep=""), nrow(combineddf))
  return(combineddf)
}, mc.cores=ncores)
allclusterdf<-do.call(rbind, newlist)
write.table(allclusterdf, file=paste(p,"average_para.tab",sep="" ), col.names = T, row.names = F, sep = "\t", quote=F)

#allclusterdf<-read.csv("/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/matrixfiles/average_para.tab", sep="\t", header=T)

##Plotting
colors=c("#0531FF","#8299FF","#000000")# wt, lean, giant colors

if (rp=="TSS"){
  p<-ggplot(allclusterdf, aes(x=BinNum, y=AvgValue))+
    geom_line(aes(color=Type))+
    theme_minimal()+
    facet_wrap(~cluster)+
    theme(axis.title.x=element_text(face="bold", size = 18),axis.title.y=element_text(face="bold", size = 18))+
    labs(x="Gene Distance, kb", y="Average Enrichment", title = paste(histonemark," ",rt2, " Regions Profile", sep=""))+
    theme(legend.text = element_text(size=14), legend.title = element_text(size = 16, face = "bold"), plot.title = element_text(hjust=0.5), strip.text.x = element_text(size=14, face="bold"))+
    scale_color_manual(values=colors, name='Phenotype')+
    geom_vline(xintercept = 0, colour="grey")+
    geom_text(aes(x=0, label="\nTSS", y=0), colour="grey", angle=90, size=2)
}

if (rp=="center"){
  p<-ggplot(allclusterdf, aes(x=BinNum, y=AvgValue))+
    geom_line(aes(color=Type))+
    theme_minimal()+
    facet_wrap(~cluster)+
    theme(axis.title.x=element_text(face="bold", size = 18),axis.title.y=element_text(face="bold", size = 18))+
    labs(x="Gene Distance, kb", y="Average Enrichment", title =paste(histonemark," ",rt2, " Regions Profile", sep=""))+
    theme(legend.text = element_text(size=14), legend.title = element_text(size = 16, face = "bold"), plot.title = element_text(hjust=0.5), strip.text.x = element_text(size=14, face="bold"))+
    scale_color_manual(values=colors, name='Phenotype')+
    geom_vline(xintercept = 0, colour="grey")+
    geom_text(aes(x=0, label="\nPC", y=0), colour="grey", angle=90, size=2)
}
ggsave(filename=paste(p3,"averageplot_",histonemark,rp,"_",rt1,"_","_para",plot_ext, sep=""), units="in", width=8, height=6)

#promoter
##Average plot
#This script will create a matrix and a line plot that takes the average values of a specific phenotype and will plot them based on the observed differentially regulated clusters
rm(list=ls())
#functions
bednames<-function(bedmat){#this function takes a BED3 file and returns a single character vector of names for regions
  bedmat[,1]=paste(bedmat[,1],rep(":", nrow(bedmat)),sep="")
  bedmat[,2]=paste(as.character(bedmat[,2]), as.character(bedmat[,3]),sep="-")
  bednames<-paste(bedmat[,1],bedmat[,2],sep="")
}

rm(list=ls())

bednames<-function(bedmat){#this function takes a BED3 file and returns a single character vector of names for regions
  bedmat[,1]=paste(bedmat[,1],rep(":", nrow(bedmat)),sep="")
  bedmat[,2]=paste(as.character(bedmat[,2]), as.character(bedmat[,3]),sep="-")
  bednames<-paste(bedmat[,1],bedmat[,2],sep="")
}
#libraries
library(ggplot2)
library(reshape2)
library(parallel)

#parameters
ncores=28
before=2000
after=1000
histonemark="H3K27ac"
plot_ext=".png"
rt1="promoter"#region type.
machine="server"
dt="seqdepthnorm"

if (rt1=="promoter"){
  rp="TSS"
  rt2="Promoter"
}
if (rt1=="enhancer"){
  rp="center"#reference point
  rt2="Enhancer"
}
#paths
if (machine=="local"){
  p="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/matrixfiles/"
  p2="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/bedfiles/"
  p3="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/plots/"
  p4="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/mousedata.tab"
}
if (machine=="server"){
  p=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",histonemark,"_",dt,"/matrixfiles/",sep="")
  p2=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",histonemark,"_",dt,"/bedfiles/",sep="")
  p3=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",histonemark,"_",dt,"/plots/",sep="")
  p4="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/mousedata.tab"
}
#regions and cluster numbers, feature data
features<-read.table(p4, header=T, sep="\t")
bedfiles<-list.files(path=p2, pattern=rt1)

if (rt1=="enhancer"){
  allbedfiles<-mclapply(1:length(bedfiles),function(i){
    test<-cbind(bednames(read.csv(paste(p2,bedfiles[i],sep=""), header = F, sep="\t")), rep(strsplit(strsplit(bedfiles[i], split="r")[[1]][3], split=".b")[[1]][1]))
  }, mc.cores = ncores)
}
if (rt1=="promoter"){
  allbedfiles<-mclapply(1:length(bedfiles),function(i){
    test<-cbind(bednames(read.csv(paste(p2,bedfiles[i],sep=""), header = F, sep="\t")), rep(strsplit(strsplit(bedfiles[i], split="r")[[1]][4], split=".b")[[1]][1]))
  }, mc.cores = ncores)
}

clusterdf<-as.data.frame(do.call(rbind, allbedfiles))
colnames(clusterdf)<-c("names", "cluster")
clusterdf<-clusterdf[which(clusterdf$names%in%as.character(read.table(paste(p2,rp,"_commonbed.bed",sep=""))$x)),]

ptvec=c("Wildtype", "Lean Nnat", "Huge Nnat")#, "Huge NnatD9")
clustervec<-levels(factor(clusterdf$cluster))

newlist<-mclapply(1:length(clustervec), function(l){
  cluster=clustervec[l]
  newclusterdf<-clusterdf[which(clusterdf$cluster==cluster),]
  gene_num=as.character(nrow(newclusterdf))

  #matrices and bedfile
  matfiles=list.files(path = p, pattern=rp)
  matfiles=matfiles[grep(pattern=".tab", x=matfiles)]
  matfiles=matfiles[grep(pattern="common", x=matfiles)]

  combined<-mclapply(1:length(matfiles), function(i){
    filen=matfiles[i]
    n=paste(strsplit(filen, split="_")[[1]][4])#makes name 
    tb=read.table(paste(p,filen,sep=""),sep="\t",stringsAsFactors = FALSE)#reads matrix file
    tb<-tb[which(rownames(tb)%in%clusterdf$names[which(clusterdf$cluster==cluster)]),]
  }, mc.cores = ncores/2)
  names<-unlist(mclapply(1:length(matfiles), function(i){
    filen=matfiles[i]
    n=paste(strsplit(filen, split="_")[[1]][4])#makes name 
  }, mc.cores = ncores/2))
  names(combined)<-names
  combinedlist<-mclapply(1:length(ptvec), function(j){
    pt=ptvec[j]
    typenums<-features$mousenumber[which(features$genopheno==pt)]
    type<-combined[which(names(combined)%in%as.character(typenums))]#returns a list of dataframes that contains cluster regions and binnned values
    typedf<-do.call(rbind,type)#combines list into a matrix
    typeavglist<-mclapply(1:length(newclusterdf$names), function(i){
      matchname<-as.character(newclusterdf$names[i])
      grep(matchname, rownames(typedf))
      filtdf<-colMeans(typedf[grep(matchname, rownames(typedf)),])
      }, mc.cores = ncores/4)
    names(typeavglist)<-as.character(newclusterdf$names)
    typeavgdf<-do.call(rbind,typeavglist)#dataframe that contains an average region enrichment value for all replicates per cluster
    typeavgvec<-colMeans(typeavgdf, na.rm = T)
    df<-rbind(rep(pt, ncol(typeavgdf)),typeavgvec, seq(-before/1000, after/1000, length.out=length(typeavgvec)))
    row.names(df)<-c("Type", "AvgValue", "BinNum")
    return(t(df))
  }, mc.cores=ncores/2)
  combineddf<-as.data.frame(do.call(rbind, combinedlist))
  combineddf$AvgValue<-as.numeric(as.character(combineddf$AvgValue))
  combineddf$BinNum<-as.numeric(as.character(combineddf$BinNum))
  combineddf$cluster<-rep(paste("Cluster ",cluster,", [",gene_num,"]",sep=""), nrow(combineddf))
  return(combineddf)
}, mc.cores=ncores)
allclusterdf<-do.call(rbind, newlist)
write.table(allclusterdf, file=paste(p,"average_para.tab",sep="" ), col.names = T, row.names = F, sep = "\t", quote=F)

#allclusterdf<-read.csv("/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/matrixfiles/average_para.tab", sep="\t", header=T)

##Plotting
colors=c("#0531FF","#8299FF","#000000")# wt, lean, giant colors

if (rp=="TSS"){
  p<-ggplot(allclusterdf, aes(x=BinNum, y=AvgValue))+
    geom_line(aes(color=Type))+
    theme_minimal()+
    facet_wrap(~cluster)+
    theme(axis.title.x=element_text(face="bold", size = 18),axis.title.y=element_text(face="bold", size = 18))+
    labs(x="Gene Distance, kb", y="Average Enrichment", title = paste(histonemark," ",rt2, " Regions Profile", sep=""))+
    theme(legend.text = element_text(size=14), legend.title = element_text(size = 16, face = "bold"), plot.title = element_text(hjust=0.5), strip.text.x = element_text(size=14, face="bold"))+
    scale_color_manual(values=colors, name='Phenotype')+
    geom_vline(xintercept = 0, colour="grey")+
    geom_text(aes(x=0, label="\nTSS", y=0), colour="grey", angle=90, size=2)
}

if (rp=="center"){
  p<-ggplot(allclusterdf, aes(x=BinNum, y=AvgValue))+
    geom_line(aes(color=Type))+
    theme_minimal()+
    facet_wrap(~cluster)+
    theme(axis.title.x=element_text(face="bold", size = 18),axis.title.y=element_text(face="bold", size = 18))+
    labs(x="Gene Distance, kb", y="Average Enrichment", title =paste(histonemark," ",rt2, " Regions Profile", sep=""))+
    theme(legend.text = element_text(size=14), legend.title = element_text(size = 16, face = "bold"), plot.title = element_text(hjust=0.5), strip.text.x = element_text(size=14, face="bold"))+
    scale_color_manual(values=colors, name='Phenotype')+
    geom_vline(xintercept = 0, colour="grey")+
    geom_text(aes(x=0, label="\nPC", y=0), colour="grey", angle=90, size=2)
}
ggsave(filename=paste(p3,"averageplot_",histonemark,rp,"_",rt1,"_","_para",plot_ext, sep=""), units="in", width=8, height=6)
