import os

counter=0

rt="promoter"#defines region type
hm="H3K27ac"#defines histone mark
bwtype=".filtered.seq_depth_norm"
dt="seqdepthnorm"
#scorefiles directory
sd="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/DNAmapping2/output/bamCoverage/"
fq1=[r for r in os.listdir(sd) if hm+bwtype in r]

af='1000'
be='2000'

#defines output directory for matrix files
om="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/matrixfiles/"

#defines output script directory
sh_p="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/scripts/"

#defines output heatmap directory
ohm="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/plots/"

#input bedfiles directory
bd0="/home/nick.wawee/nickS/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/bedfiles/"
bd=bd0+rt+"_*"
#output bedfiles directory
bdo="/home/nick.wawee/nickS/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/bedfiles/"

#reference point, lower case c for center
rp='TSS'
#reference point label
rpl='TSS'

#binsize
bs='10'

#number of cores
nc='4'

#regions labeling (length should equal to the number of input bedfiles)
#rl=[]
#beds=[bed for bed in os.listdir(bd0) if 'enhancer' in bed]
#for bed in beds:
#    rl.append(bed.split('r')[2].split(".b")[0])
#rl='" "'.join(rl)
rl=' "1" "2" "3" "4" "5" "6" "7" "9"'
for name in fq1:
    if "lean_nnat" in name:
        f=open(sh_p+rt+str(counter)+".sh","w")
        f.write('computeMatrix reference-point -R '+bd+' -S '+sd+name+' -b '+be+' -a '+af+' --referencePoint '+rp+' --missingDataAsZero --skipZeros -bs '+bs+' -o '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz  --outFileNameMatrix '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.tab -p '+nc+'\nplotHeatmap -m '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz --legendLocation none --colorList "white,blue" --alpha 0.5  --zMin 0 --zMax 4 --refPointLabel '+rpl+' --interpolationMethod bilinear --boxAroundHeatmaps no --xAxisLabel " " --samplesLabel M'+name.split(sep="_")[0]+ ' --whatToShow "heatmap only" --sortRegions no --regionsLabel'+rl+' -o '+ohm+'heatmap_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.png --outFileSortedRegions '+bdo+'sortedregions_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.bed')
        f.close()
    if "huge" in name:
        f=open(sh_p+rt+str(counter)+".sh","w")
        f.write('computeMatrix reference-point -R '+bd+' -S '+sd+name+' -b '+be+' -a '+af+' --referencePoint '+rp+' --missingDataAsZero --skipZeros -bs '+bs+' -o '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz  --outFileNameMatrix '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.tab -p '+nc+'\nplotHeatmap -m '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz --legendLocation none --colorList "white,blue" --zMin 0 --zMax 4 --refPointLabel ' +rpl+' --interpolationMethod bilinear --boxAroundHeatmaps no --xAxisLabel " " --samplesLabel M'+name.split(sep="_")[0]+ ' --whatToShow "heatmap only" --sortRegions no --refPointLabel '+rpl+' --regionsLabel'+rl+' -o '+ohm+'heatmap_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.png --outFileSortedRegions '+bdo+'sortedregions_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.bed')
        f.close()
    if "lean_wt" in name:
        f=open(sh_p+rt+str(counter)+".sh","w")
        f.write('computeMatrix reference-point -R '+bd+' -S '+sd+name+' -b '+be+' -a '+af+' --referencePoint '+rp+' --missingDataAsZero --skipZeros -bs '+bs+' -o '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz  --outFileNameMatrix '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.tab -p '+nc+'\nplotHeatmap -m '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz --legendLocation none --colorList "white,black" --zMin 0 --zMax 4 --interpolationMethod bilinear --boxAroundHeatmaps no --xAxisLabel " " --samplesLabel M'+name.split(sep="_")[0]+ ' --whatToShow "heatmap only" --sortRegions no --refPointLabel '+rpl+' --regionsLabel'+rl+' -o '+ohm+'heatmap_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.png --outFileSortedRegions '+bdo+'sortedregions_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.bed')
        f.close()
    if "lean_wt" and "687" in name:
        f=open(sh_p+rt+str(counter)+".sh","w")
        f.write('computeMatrix reference-point -R '+bd+' -S '+sd+name+' -b '+be+' -a '+af+' --referencePoint '+rp+' --missingDataAsZero --skipZeros -bs '+bs+' -o '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz  --outFileNameMatrix '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.tab -p '+nc+'\nplotHeatmap -m '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz --legendLocation none --colorList "white,black" --zMin 0 --zMax 4 --interpolationMethod bilinear --boxAroundHeatmaps no --xAxisLabel " " --refPointLabel '+rpl+' --samplesLabel M'+name.split(sep="_")[0]+ ' --whatToShow "heatmap and colorbar" --sortRegions no  --regionsLabel'+rl+' -o '+ohm+'heatmap_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.png --outFileSortedRegions '+bdo+'sortedregions_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.bed')
        f.close()
    if "lean_nnat" and "686" in name:
        f=open(sh_p+rt+str(counter)+".sh","w")
        f.write('computeMatrix reference-point -R '+bd+' -S '+sd+name+' -b '+be+' -a '+af+' --referencePoint '+rp+' --missingDataAsZero --skipZeros -bs '+bs+' -o '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz  --outFileNameMatrix '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.tab -p '+nc+'\nplotHeatmap -m '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz --legendLocation none --colorList "white,blue" --alpha 0.5  --zMin 0 --zMax 4 --interpolationMethod bilinear --boxAroundHeatmaps no --xAxisLabel " " --samplesLabel M'+name.split(sep="_")[0]+ ' --whatToShow "heatmap and colorbar" --sortRegions no --refPointLabel '+rpl+' --regionsLabel'+rl+' -o '+ohm+'heatmap_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.png --outFileSortedRegions '+bdo+'sortedregions_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.bed')
        f.close()
    if "huge" and "680" in name:
        f=open(sh_p+rt+str(counter)+".sh","w")
        f.write('computeMatrix reference-point -R '+bd+' -S '+sd+name+' -b '+be+' -a '+af+' --referencePoint '+rp+' --missingDataAsZero --skipZeros -bs '+bs+' -o '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz  --outFileNameMatrix '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.tab -p '+nc+'\nplotHeatmap -m '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz --legendLocation none --colorList "white,blue" --zMin 0 --zMax 4 --interpolationMethod bilinear --boxAroundHeatmaps no --xAxisLabel " " --samplesLabel M'+name.split(sep="_")[0]+ ' --whatToShow "heatmap and colorbar" --sortRegions no --refPointLabel '+rpl+' --regionsLabel'+rl+' -o '+ohm+'heatmap_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.png --outFileSortedRegions '+bdo+'sortedregions_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.bed')
        f.close()
    counter+=1


import os
#enhancers
counter=0

rt="enhancer"#defines region type
hm="H3K27ac"#defines histone mark
bwtype=".filtered.seq_depth_norm"
dt="seqdepthnorm"
#scorefiles directory
sd="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/DNAmapping2/output/bamCoverage/"
fq1=[r for r in os.listdir(sd) if hm+bwtype in r]

af='5000'
be='5000'

#defines output directory for matrix files
om="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/matrixfiles/"

#defines output script directory
sh_p="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/scripts/"

#defines output heatmap directory
ohm="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/plots/"

#input bedfiles directory
bd0="/home/nick.wawee/nickS/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/bedfiles/"
bd=bd0+rt+"_*"
#output bedfiles directory
bdo="/home/nick.wawee/nickS/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/bedfiles/"

#reference point, lower case c for center
rp='center'
#reference point label
rpl='Center'

#binsize
bs='10'

#number of cores
nc='4'

#regions labeling (length should equal to the number of input bedfiles)
#rl=[]
#beds=[bed for bed in os.listdir(bd0) if 'enhancer' in bed]
#for bed in beds:
#    rl.append(bed.split('r')[2].split(".b")[0])
#rl='" "'.join(rl)
rl=' "1" "2" "3" "4" "5" "6" "8" "9"'
for name in fq1:
    if "lean_nnat" in name:
        f=open(sh_p+rt+str(counter)+".sh","w")
        f.write('computeMatrix reference-point -R '+bd+' -S '+sd+name+' -b '+be+' -a '+af+' --referencePoint '+rp+' --missingDataAsZero --skipZeros -bs '+bs+' -o '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz  --outFileNameMatrix '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.tab -p '+nc+'\nplotHeatmap -m '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz --legendLocation none --colorList "white,blue" --alpha 0.5  --zMin 0 --zMax 4 --refPointLabel '+rpl+' --interpolationMethod bilinear --boxAroundHeatmaps no --xAxisLabel " " --samplesLabel M'+name.split(sep="_")[0]+ ' --whatToShow "heatmap only" --sortRegions no --regionsLabel'+rl+' -o '+ohm+'heatmap_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.png --outFileSortedRegions '+bdo+'sortedregions_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.bed')
        f.close()
    if "huge" in name:
        f=open(sh_p+rt+str(counter)+".sh","w")
        f.write('computeMatrix reference-point -R '+bd+' -S '+sd+name+' -b '+be+' -a '+af+' --referencePoint '+rp+' --missingDataAsZero --skipZeros -bs '+bs+' -o '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz  --outFileNameMatrix '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.tab -p '+nc+'\nplotHeatmap -m '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz --legendLocation none --colorList "white,blue" --zMin 0 --zMax 4 --refPointLabel ' +rpl+' --interpolationMethod bilinear --boxAroundHeatmaps no --xAxisLabel " " --samplesLabel M'+name.split(sep="_")[0]+ ' --whatToShow "heatmap only" --sortRegions no --refPointLabel '+rpl+' --regionsLabel'+rl+' -o '+ohm+'heatmap_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.png --outFileSortedRegions '+bdo+'sortedregions_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.bed')
        f.close()
    if "lean_wt" in name:
        f=open(sh_p+rt+str(counter)+".sh","w")
        f.write('computeMatrix reference-point -R '+bd+' -S '+sd+name+' -b '+be+' -a '+af+' --referencePoint '+rp+' --missingDataAsZero --skipZeros -bs '+bs+' -o '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz  --outFileNameMatrix '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.tab -p '+nc+'\nplotHeatmap -m '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz --legendLocation none --colorList "white,black" --zMin 0 --zMax 4 --interpolationMethod bilinear --boxAroundHeatmaps no --xAxisLabel " " --samplesLabel M'+name.split(sep="_")[0]+ ' --whatToShow "heatmap only" --sortRegions no --refPointLabel '+rpl+' --regionsLabel'+rl+' -o '+ohm+'heatmap_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.png --outFileSortedRegions '+bdo+'sortedregions_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.bed')
        f.close()
    if "lean_wt" and "687" in name:
        f=open(sh_p+rt+str(counter)+".sh","w")
        f.write('computeMatrix reference-point -R '+bd+' -S '+sd+name+' -b '+be+' -a '+af+' --referencePoint '+rp+' --missingDataAsZero --skipZeros -bs '+bs+' -o '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz  --outFileNameMatrix '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.tab -p '+nc+'\nplotHeatmap -m '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz --legendLocation none --colorList "white,black" --zMin 0 --zMax 4 --interpolationMethod bilinear --boxAroundHeatmaps no --xAxisLabel " " --refPointLabel '+rpl+' --samplesLabel M'+name.split(sep="_")[0]+ ' --whatToShow "heatmap and colorbar" --sortRegions no  --regionsLabel'+rl+' -o '+ohm+'heatmap_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.png --outFileSortedRegions '+bdo+'sortedregions_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.bed')
        f.close()
    if "lean_nnat" and "686" in name:
        f=open(sh_p+rt+str(counter)+".sh","w")
        f.write('computeMatrix reference-point -R '+bd+' -S '+sd+name+' -b '+be+' -a '+af+' --referencePoint '+rp+' --missingDataAsZero --skipZeros -bs '+bs+' -o '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz  --outFileNameMatrix '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.tab -p '+nc+'\nplotHeatmap -m '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz --legendLocation none --colorList "white,blue" --alpha 0.5  --zMin 0 --zMax 4 --interpolationMethod bilinear --boxAroundHeatmaps no --xAxisLabel " " --samplesLabel M'+name.split(sep="_")[0]+ ' --whatToShow "heatmap and colorbar" --sortRegions no --refPointLabel '+rpl+' --regionsLabel'+rl+' -o '+ohm+'heatmap_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.png --outFileSortedRegions '+bdo+'sortedregions_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.bed')
        f.close()
    if "huge" and "680" in name:
        f=open(sh_p+rt+str(counter)+".sh","w")
        f.write('computeMatrix reference-point -R '+bd+' -S '+sd+name+' -b '+be+' -a '+af+' --referencePoint '+rp+' --missingDataAsZero --skipZeros -bs '+bs+' -o '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz  --outFileNameMatrix '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.tab -p '+nc+'\nplotHeatmap -m '+om+'refpoint_'+rp+'_'+name.split(sep=".")[0]+'_p'+be+'m'+af+'_'+'.gz --legendLocation none --colorList "white,blue" --zMin 0 --zMax 4 --interpolationMethod bilinear --boxAroundHeatmaps no --xAxisLabel " " --samplesLabel M'+name.split(sep="_")[0]+ ' --whatToShow "heatmap and colorbar" --sortRegions no --refPointLabel '+rpl+' --regionsLabel'+rl+' -o '+ohm+'heatmap_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.png --outFileSortedRegions '+bdo+'sortedregions_refpoint_'+rp+'_p'+be+'m'+af+'_'+name.split(sep=".")[0]+'.bed')
        f.close()
    counter+=1

