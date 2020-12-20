import os

hm="H3K27ac"
dt="seqdepthnorm"
counter=0
#retrieves all input file names, bigwig
p1="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/DNAmapping2/output/bamCoverage/"
fq1=[r for r in os.listdir(p1) if hm+".filtered.seq_depth_norm" in r]

#defines output directory
op="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/bedfiles/"

#defines output script directory
sh_p="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/scripts/"

#Rscriptpath
rsp="/primary/projects/mnp/tools/anaconda3/bin/Rscript"

#Running R script path
rrsp="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/Get.enrichment.forUnions.r"

#bed file path
bedpath="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/Peak_Counting_Compare/"+hm+"_merged.bed"

#number of cores
ncores="20"

for name in fq1:
    print(rsp+ " "+ rrsp+ " " +bedpath + " "+ p1+name +" "+op+ " "+ ncores)
    f=open(sh_p+"e"+str(counter)+".sh","w")
    f.write(rsp+ " "+ rrsp+ " " +bedpath + " "+ p1+name +" "+op+name.split(sep=".bw")[0]+".bed"+ " "+ ncores)
    f.close()
    counter+=1
