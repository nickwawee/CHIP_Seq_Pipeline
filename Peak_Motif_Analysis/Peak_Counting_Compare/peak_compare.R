##Finding unique peaks
#This script will find unique peaks with each phenotype as well as between phenotypes
rm(list=ls())
library(GenomicRanges)
library(ChIPpeakAnno)
library(eulerr)
library(limma)
library(argparser, quietly=TRUE)#use when you dont have write permissions, lib.loc="/primary/projects/mnp/nick/software/R/lib")

# Create a parser
p <- arg_parser("Create a venn diagram that shows how different peaks within and between phenotypes.
                , it assumes that the file names are in the following format: animalnumber_phenotype_genotype_sex_antibody. 
                It also assumes that there are six replicates for the name lean_wt and four replicates for the names huge_nnat and lean_nnat.
                It filters out JH and BL chromosomes. It also assumes three phenotype_genotype combinations.")

# Add command line arguments
p <- add_argument(p, "--bp", help="bed path, path to bed files", type="character")
p <- add_argument(p, "--pp", help="plot path, path where plot will save to. Must have / at end", type="character")
p <- add_argument(p, "--be", help="bed extension, bed file extension produced by MACS2. ex.) .filtered.BAMPE_peaks.broadPeak ", type="character")
p <- add_argument(p, "--p", help="number of processors", type="integer")
p <- add_argument(p, "--mo", help="minimum overlap, number of basepairs that can overlap to be considered at peak", type="integer")
# Parse the command line arguments
argv <- parse_args(p)

bp=argv$bp
od=argv$pp
bedext=argv$be
nc=argv$p
mo=argv$mo

bedfiles=list.files(path=bp, pattern=bedext)

#finding all antibodies
splitby=paste(strsplit2(bedext,split="")[,1:2], collapse="")
hms=strsplit2(levels(factor(strsplit2(bedfiles,"_")[,5])),split=splitby)[,1]

#finding all phenotypes
ptvec=levels(factor(paste(strsplit2(bedfiles,"_")[,2], strsplit2(bedfiles,"_")[,3],sep = "_")))

mclapply(1:length(hms), function(k){
  hm=hms[k]
  grlist_pt=mclapply(1:length(ptvec), function(i){#this will return a list of GRanges objects that contains each pt_genotype combination that is common between each replicate
    pt=ptvec[i]
    bfs_p=bedfiles[grep(pattern=pt, x=bedfiles)]
    grlist=mclapply(1:length(bfs_p), function(j){#this will return a granges object for each replicate
      bed=read.table(file=paste(bp,bfs_p[j], sep=""), sep="\t")
      if (length(grep('JH', bed$V1))>0){bed=bed[-grep('JH', bed$V1),]}
      if (length(grep('GL', bed$V1))>0){bed=bed[-grep('GL', bed$V1),]}
      GRanges(seqnames = as.character(bed$V1), ranges=IRanges(start=as.integer(bed$V2), end=as.integer(bed$V3), names=as.character(bed$V4)))
    }, mc.cores=nc/4)
    if (pt=="lean_wt"){
      ol=findOverlapsOfPeaks(grlist[[1]], grlist[[2]], grlist[[3]], grlist[[4]], grlist[[5]], minoverlap = mo)
      ol2=findOverlapsOfPeaks(ol[["peaklist"]][["grlist..1..///grlist..2..///grlist..3..///grlist..4..///grlist..5.."]], grlist[[6]], minoverlap = mo)
      return(ol2[["peaklist"]][["ol...peaklist......grlist..1.....grlist..2.....grlist..3.....grlist..4.....grlist..5.....///grlist..6.."]])
    }
    else if (pt=="lean_nnat" | pt=="huge_nnat"){
      ol=findOverlapsOfPeaks(grlist[[1]], grlist[[2]], grlist[[3]], grlist[[4]], minoverlap = mo)
      return(ol[["peaklist"]][["grlist..1..///grlist..2..///grlist..3..///grlist..4.."]])
    }
    else{stop("Names do not have lean_wt, huge_nnat, or lean_nnat in them. Please change code and appropriate replicate number accordingly.")}
  }, mc.cores = nc/2)

  names(grlist_pt)<-ptvec
  save(grlist_pt, file=paste(od,hm,"_phenopeaks.rda", sep=""))
  #next, the interphenotype peaks will be compared and plotted via a venn diagram
  ol3=findOverlapsOfPeaks(grlist_pt[[1]], grlist_pt[[2]], grlist_pt[[3]], minoverlap = mo)

  #getting numbers of peaks
  h_nnat=as.numeric(length(ol3[["peaklist"]][["grlist_pt..3.."]]@ranges))
  l_nnat=as.numeric(length(ol3[["peaklist"]][["grlist_pt..2.."]]@ranges))
  h_l_nnat=as.numeric(length(ol3[["peaklist"]][["grlist_pt..2..///grlist_pt..3.."]]@ranges))
  l_wt=as.numeric(length(ol3[["peaklist"]][["grlist_pt..1.."]]@ranges))
  l_wt_h_nnat=as.numeric(length(ol3[["peaklist"]][["grlist_pt..1..///grlist_pt..3.."]]@ranges))
  l_wt_l_nnat=as.numeric(length(ol3[["peaklist"]][["grlist_pt..1..///grlist_pt..2.."]]@ranges))
  l_wt_h_l_nnat=as.numeric(length(ol3[["peaklist"]][["grlist_pt..1..///grlist_pt..2..///grlist_pt..3.."]]@ranges))

  colors=c( "#0531FF", "#8299FF","#000000")# wt, lean, giant colors
  euobj=euler(c(`Huge Nnat`=h_nnat, `Lean Nnat`=l_nnat,"Huge Nnat&Lean Nnat"=h_l_nnat, `WT`=l_wt, "WT&Huge Nnat"=l_wt_h_nnat, "WT&Lean Nnat"=l_wt_l_nnat, "WT&Lean Nnat&Huge Nnat"= l_wt_h_l_nnat))
  png(filename=paste(od,hm,"_interphenodiagram.png", sep=""), units='in', height=6, width=6, res=600)
  plot(euobj,quantities=list(fontsize=12, col="white"), col=colors, fill=colors,cat=list(col=rep('white',3)), labels=list(col=rep('white',3)))
  dev.off()
}, mc.cores = nc)
