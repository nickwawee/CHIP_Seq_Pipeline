##peak counting
#This script will count peaks from bed files and plot them
rm(list=ls())
library(parallel)
library(ggplot2)
library(Hmisc)
library(argparser, quietly=TRUE)#use when you dont have write permissions, lib.loc="/primary/projects/mnp/nick/software/R/lib")

# Create a parser
p <- arg_parser("Create a plot that counts peaks, it assumes that the file names are in the following format: animalnumber_phenotype_genotype_sex_antibody. It filters out JH and BL chromosomes")

# Add command line arguments
p <- add_argument(p, "--bp", help="bed path, path to bed files", type="character")
p <- add_argument(p, "--pp", help="plot path, path where plot will save to. Must have / at end", type="character")
p <- add_argument(p, "--be", help="bed extension, bed file extension produced by MACS2. ex.) .filtered.BAMPE_peaks.broadPeak ", type="character")
p <- add_argument(p, "--p", help="number of processors", type="integer")
# Parse the command line arguments
argv <- parse_args(p)

bp=argv$bp
plotp=argv$pp
bedext=argv$be
nc=argv$p

bedfiles=list.files(path=bp, pattern=bedext)
counts=mclapply(1:length(bedfiles), function(i){
  bedfile=read.table(file=paste(bp,bedfiles[i], sep=""), sep="\t")
  if (length(grep('JH', bedfile$V1)>0)){
    bedfile=bedfile[-grep('JH', bedfile$V1),]
  }
  if (length(grep('GL', bedfile$V1))){
    bedfile=bedfile[-grep('GL', bedfile$V1),]
  }
  
  return(data.frame(name=as.character(bedfiles[i]), pc=as.integer(nrow(bedfile)), gt=paste(unlist(strsplit(bedfiles[i], split="_"))[2:3], collapse = "_"), pt=unlist(strsplit(bedfiles[i], split="_"))[2], hm=unlist(strsplit(unlist(strsplit(bedfiles[i], split=".f"))[1], split="_"))[5]))
}, mc.cores=nc)

countdf<-do.call(rbind,counts)
if (length(grep('Error', x=countdf$pc))>0){
   countdf<-countdf[-grep('Error', x=countdf$pc),]#removes 0 counts
}
countdf$pc<-as.numeric(countdf$pc)

##plotting

p2=ggplot(countdf, aes(x=gt, y=pc))+
  geom_dotplot(binaxis = "y", stackdir="center")+ 
  theme_minimal()+
  facet_wrap(~hm, scales="free_y")+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="red", width=0.2) +
  stat_summary(fun.y=mean, geom="point", color="red")+
  labs(x='Phenotype_Genotype', y="# of Peaks", title='MACS2 Peaks Called')+
  theme(strip.text = element_text(face='bold', size=12), axis.title=element_text(face='bold', size=12), plot.title = element_text(hjust=0.5, face='bold', size=16), axis.text.x=element_text(angle=25))

ggsave(filename=paste(plotp,"PeakCount.png",sep=""), dpi=600, plot=p2, height=5, width=7, units="in")
