##PCA(s) on peak AUC's
#This script will perform and plot PCA results from the Peak AUCS
rm(list=ls())
library(ggplot2)
library(ggpubr)

percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
histonemark="H3K27ac"
dt="seqdepthnorm"
machine="server"
if (machine=="local"){
   load(paste("/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",histonemark,"_",dt,"/matrixfiles/enrichmentdf.RData", sep=""))
   features=read.table(file=paste("/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/mousedata.tab", sep=""), sep="\t", header=T)
}else{
   load(paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",histonemark,"_",dt,"/matrixfiles/enrichmentdf.RData", sep=""))
   features=read.table(file=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/mousedata.tab", sep=""), sep="\t", header=T)
}
enrichmat<-enrichmentdf[,-ncol(enrichmentdf)]
features<-features[order(features$mousenumber),]
filter='no'
badind<-as.numeric()
d9='no'

scaledmat<-scale(enrichmat, center=TRUE, scale=TRUE)

#no d9
if (d9=='no'){
  if (length(grep('d9',colnames(scaledmat)))>0){scaledmat=scaledmat[,-grep('d9',colnames(scaledmat))]}
  features=features[-which(features$genopheno=='Huge NnatD9'),]
}
pcares<-prcomp(t(scaledmat), center=TRUE, scale.=FALSE)
eigvals=pcares$sdev^2#calculates eigenvalues from the PCA
pc1_var=eigvals[1]/sum(eigvals)#calculates percent variance for PC1
pc2_var=eigvals[2]/sum(eigvals)
pc3_var=eigvals[3]/sum(eigvals)
pc1_percvar=percent(pc1_var)
pc2_percvar=percent(pc2_var)

df=data.frame(PC1=pcares$x[,1], PC2=pcares$x[,2])
df<-cbind(df,features)

colors=c( "#0531FF", "#8299FF","#000000")# three group color palette

p1=ggplot(df,aes(x=PC1, y=PC2))+
  geom_point(aes(color=pheno))+
  theme_minimal()+
  theme(axis.title.x=element_text(face="bold", size = 18),axis.title.y=element_text(face="bold", size = 18))+
  labs(x=paste("PC1", pc1_percvar), y=paste("PC2",pc2_percvar))+
  theme(legend.text = element_text(size=16), legend.title = element_text(size = 18, face = "bold"))+
  theme(legend.position = c(1, 0),legend.justification = c(1, -.8), legend.background = element_rect(colour = NA, fill = "white"))+
  scale_color_manual(values=colors, name='Phenotype')+ 
  ggtitle('')

p2=ggplot(df,aes(x=pheno, y=PC1))+
  geom_dotplot(binaxis = "y", stackdir="center")+ 
  theme_minimal()+
  theme(axis.title.x=element_text(face="bold", size = 18),axis.title.y=element_text(face="bold", size = 18))+
  labs(x="Phenotype", y=paste("PC1",pc1_percvar))+
  scale_color_discrete(name='Genotype Phenotype')+
  theme(legend.text = element_text(size=16), legend.title = element_text(size = 18, face = "bold"))+
  theme(legend.position = "none",legend.justification = c(1, -2), legend.background = element_rect(colour = NA, fill = "white"))+
  scale_color_manual(values=colors)+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="red", width=0.2) +
  stat_summary(fun.y=mean, geom="point", color="red")

p3=ggarrange(p1,p2,nrow=1)

ggsave(p3, file=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",histonemark,"_",dt,"/plots/pca.png", sep=""), units="in", height=4, width=8, dpi=600)
