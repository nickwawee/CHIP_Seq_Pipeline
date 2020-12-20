##Convert homer results into table for GO analysis
rm(list=ls())
library(limma)
library(ggplot2)
p="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/Motif_Analysis_4/homer/knownResults.txt"
loaddf<-read.csv(file=p, sep="\t", header = T)

goodpvals<-loaddf$P.value[grep(')/', loaddf$Motif.Name)]
goodnames<-as.character(loaddf$Motif.Name[grep(')/', loaddf$Motif.Name)])
genenames<-strsplit2(goodnames, split="[()]")[,1]

df<-cbind(genenames,goodpvals)[-which(goodpvals>0.05),]
write.table(df, file="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/Motif_Analysis/homer/diff_giantvswt/pvaltable.txt",col.names = F, row.names = F, quote = F)



gotermtab<-read.csv(file="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/Motif_Analysis/homer/diff_giantvswt/go_res.csv")
tidygo<-gotermtab[-which(gotermtab$adjustedPValue_giantvswt_homer>0.05),]
length(grep("pancrea", tidygo$description_giantvswt_homer))
length(grep("lipid", tidygo$description_giantvswt_homer))

#making volcano plot
volcanodf<-data.frame(Log2fc=log2(as.numeric(sub("%","",loaddf$X..of.Target.Sequences.with.Motif))/as.numeric(sub("%","",loaddf$X..of.Background.Sequences.with.Motif))), neglogp=-1*log10(loaddf$P.value))

volcanodf$genes<-NA
volcanodf$threshold3=volcanodf$logFC>4 & volcanodf$pvals<pvalthresh | volcanodf$logFC< -6.5&volcanodf$pvals<pvalthresh
volcanodf$genes[which(volcanodf$threshold3)]<-tssgenes[which(volcanodf$threshold3)]

ggplot(volcanodf, aes(x=2^Log2fc, y=neglogp))+
  geom_point(aes( alpha=0.1))+
  theme_minimal()+
  labs(x=bquote(""~log[2]~"["~over(WT,Huge)~"]"), y=bquote("-"~log[10]~"[p-value]"))+
  scale_x_continuous(limits=c(-12,12))+
  theme(legend.position = "none")#+
  #geom_text_repel(aes( label=genes))