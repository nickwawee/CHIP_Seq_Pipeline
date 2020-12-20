##Make bar plot for GSEA terms
#This script will make bar plot(s) for GSEA terms
rm(list=ls())
#library(data.table)
library(parallel)
library(ggplot2)
library(limma)
library(ggrepel)

hm="H3K27ac"
dt="seqdepthnorm"
phenotype="wt"
rt="promoter"#region type
rt2="Promoter"
p=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",hm,"_",dt,"/GSEA/", sep="")
allfiles<-list.files(paste(p, sep=""))
filtfiles<-allfiles[grep(pattern="Gsea",x=allfiles)]#filters folders outputted by Gsea
filtfiles=filtfiles[grep(pattern=rt, x=filtfiles)]
#loading
dfs<-mclapply(1:length(filtfiles), function(i){
  filtfiles2<-list.files(paste(p, filtfiles[i], sep=""))
  if (rt=="promoter"){
    num=strsplit(filtfiles[i], split="a.")[[1]][2]
  }else{num=strsplit(filtfiles[i], split="a.")[[1]][3]}
  filtfiles3=filtfiles2[grep(pattern=num, x=filtfiles2)]
  xlfile=filtfiles3[grep(pattern=paste("gsea_report_for_",phenotype,"_",num,".xls",sep=""), x=filtfiles3)]
  xltabl=fread(paste(p,filtfiles[i],"/",xlfile,sep=""), drop=3)
  data=as.data.frame(cbind(as.character(xltabl$NAME),as.numeric(xltabl$`FDR q-val`),as.numeric(xltabl$`NOM p-val`), as.numeric(xltabl$NES),as.character(rep(strsplit(filtfiles[i], split="_")[[1]][2]))))
})
combined<-do.call(rbind, dfs)
names(combined)<-c("Name", "FDR", "Nom_P","NES", "Dataset")

combined$Name<-as.character(combined$Name)
combined$FDR<-as.numeric(as.character(combined$FDR))
combined$Nom_P<-as.numeric(as.character(combined$Nom_P))
combined$NES<-as.numeric(as.character(combined$NES))
combined$Dataset<-as.character(combined$Dataset)

#renaming Dataset types
combined$Dataset[which(combined$Dataset=="c2cgp")]<-"Chemical and Genetic Peturbations"
combined$Dataset[which(combined$Dataset=="c2cp")]<-"Canonical Pathways"
#combined$Dataset[which(combined$Dataset=="c5gobp")]<-"GO Biological Process"
#combined$Dataset[which(combined$Dataset=="c5gocc")]<-"GO Cellular Component"

#filtering
combined<-combined[-which(combined$FDR >0.01),]

##Plotting
#scatter plot
#p=ggplot(combined, aes(x=NES, y=-log10(FDR)))+
#  geom_point( alpha=0.1)+
#  theme_minimal()+
#  fabs(x="NES", y=bquote("-"~log[10]~"[FDR]"))+
#  scale_x_continuous(limits=c(-3,3))+
#  theme(legend.position = "none")+
#  geom_text_repel(aes( label=Name), size= 2.5)
#p

#bargraph
p1=ggplot(combined)+
  geom_col(aes(x=reorder(Name,abs(NES)), y=NES, fill=FDR))+
  facet_wrap(~Dataset,scales = "free_y")+
  #facet_grid(cols=vars(Dataset), scales = "free_y")+
  coord_flip()+
  theme_minimal()+
  theme(axis.title.y=element_blank(), axis.text = element_text(size=12), strip.text = element_text(face="bold"), plot.title=element_text(face='bold', hjust=0.5))+
  ggtitle(paste(rt2, " GSEA Results", sep=""))
p1

ggsave(filename=paste(p,rt,"_",phenotype,"_","GSEAplot.png", sep=""), device = "png",plot = p1, height =8, width=18, units = "in" )
