##GSEA file write
#This script will write expression and phenotype files required for the GSEA software
rm(list=ls())
library(limma)
library(parallel)
rt="Promoter"#region type
hm="H3K27ac"
dt="seqdepthnorm"
nc=4
#loading
load("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/H3K27ac_ChIPminusInput/matrixfiles/enrichmentdf.RData")
load("/primary/projects/mnp/tools/data/repository/organisms/GRCm38.p6_GENCODE/annotation/gene.filter/GRCm36.p6.GENCODE.RefSeq.filt.RData")

#writing expression file, .txt
prommat<-enrichmentdf[which(enrichmentdf$Type==rt), -ncol(enrichmentdf)]
bed=data.frame(chr=as.character(strsplit2(rownames(prommat), split = ":")[,1]), st=as.numeric(strsplit2(strsplit2(rownames(prommat), split = ":")[,2], split="-")[,1]), en=as.numeric(strsplit2(strsplit2(rownames(prommat), split = ":")[,2], split="-")[,2]))

if (rt=="Enhancer"){
  enhagenes<-unlist(mclapply(1:nrow(prommat), function(i){
    st=as.numeric(bed$st[i]); en=as.numeric(bed$en[i]); chr=as.character(bed$ch[i]);
    pc=ceiling((st+en)/2)#peak center
    ind=which.min(abs(pc-genes$tss))
    gene<-genes$geneName[ind]
    return(gene)
  }, mc.cores = nc))
  combined=prommat
  combined=combined[,-grep("lean_nnat", colnames(combined))]#removing lean mice
  combined<-combined[,order(strsplit2(colnames(combined),split="_")[,2])]#sorts according to genotype
  combined<-cbind(enhagenes, rep("na",nrow(combined)), combined)
  colnames(combined)<-c("NAMES", "DESCRIPTION", colnames(combined)[c(-1,-2)])
}
if (rt=="Promoter"){
    tssgenes<-mclapply(1:nrow(bed), function(i){
      st=as.numeric(bed$st[i]); en=as.numeric(bed$en[i]); chr=as.character(bed$ch[i]);
      ind=which(genes$chr==chr & genes$tss>=st & genes$tss<=en)
      bedname=paste(bed[i,1],paste(bed[i,2:3],collapse="-"), sep=":")
      df<-data.frame(genesym=as.character(genes$geneName[ind]), geneid=as.character(genes$id[ind]), bednames=as.character(rep(bedname, length(ind))))
      return(df)
    }, mc.cores = nc)
    tss_df<-do.call(rbind, tssgenes)
    duplicates<-prommat[match(tss_df[which(duplicated(tss_df$bednames)),3], rownames(prommat)),]
    row.names(duplicates)<-as.character(tss_df$genesym[which(duplicated(tss_df$bednames))])
    
    uniques<-prommat[match(tss_df[which(!duplicated(tss_df[,3])),3], rownames(prommat)),]
    row.names(uniques)<-as.character(tss_df$genesym[which(!duplicated(tss_df$bednames))])
    combined<-as.data.frame(rbind(duplicates,uniques))
    combined<-combined[,-grep("lean_nnat", colnames(combined))]#removing lean mice
    combined<-combined[,order(strsplit2(colnames(combined),split="_")[,2])]#sorts according to genotype
    combined<-cbind(rownames(combined), rep("na",nrow(combined)), combined)
    colnames(combined)<-c("NAMES", "DESCRIPTION", colnames(combined)[c(-1,-2)])
}

write.table(combined, file=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",hm,"_",dt,"/GSEA/",rt,".txt", sep=""), col.names = T, row.names = F, quote = F, sep="\t")

#writing phenotype file, .cls
npheno=2
line1=paste(c(ncol(combined)-2,npheno,1 ), collapse = " ")
line2="# huge wt"
line3=paste(strsplit2(colnames(combined), split="_")[,2][c(-1,-2)], collapse = " ")
test<-as.matrix(c(line1,line2,line3))
write.table(c(line1,line2,line3), file=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",hm,"_",dt,"/GSEA/",rt,".cls",sep=""),col.names = F, row.names = F, quote = F, sep="\t")
