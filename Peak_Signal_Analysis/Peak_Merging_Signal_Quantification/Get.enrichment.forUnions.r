library(rtracklayer)
require("parallel")

args = commandArgs(TRUE) 
if (length(args)==0){cat("please enter file names and number of cores");q("no")}
if (length(args)<4){cat("one or more argument(s) missing");q("no")}

fbd<-args[[1]] ### the input bed file  (merged peaks called using MACS2)
fbw<-args[[2]] ### the input bigwig files (generated  deeptools BamCompare directory)
fbd.out<-args[[3]] ### output file
ncores<-as.numeric(args[[4]])  ### number of cores

#reading the bed
bed=read.csv(fbd,sep="\t", header=FALSE);names(bed)=c("chr","start","end")

# read bigwig
data<-import.bw(con=fbw)
chrn<-as.vector(data@seqnames)
starts<-data@ranges@start
ends<-data@ranges@start+data@ranges@width-1
enr<-data$score

 val=unlist(mclapply(1:length(bed$chr),function(i){
        chr=bed$chr[i]; st=bed$start[i]; en=bed$end[i];
        ix=which(chrn==chr & starts>=st & starts<=en & ends>=st & ends<=en)
	r=sum(enr[ix]*data@ranges@width[ix])
	return(r)
},mc.cores=ncores))

t=cbind(as.character(bed$chr),bed$start,bed$end,val)
write.table(t,file=fbd.out,sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
cat(paste(fbd.out," saved successfully!",sep=""))
 
