##This script combines multiple heatmaps according to genotype-phenotype combination
#promoters
rm(list=ls())
library(magick)
colors=c("#000000", "#8299FF", "#0531FF")# wt, lean, giant colors
rt="Promoter"
hm="H3K27ac"
dt="seqdepthnorm"
machine="local"
if (machine=="local"){
   p=paste("/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",hm,"_",dt,"/plots/", sep="")
}else{
   p=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",hm,"_",dt,"/plots/", sep="")
}
images<-list.files(path=p)

#filtering

images<-images[grep(x=images,"TSS")]
id=paste(strsplit(images, split="_")[[1]][3:4], collapse ="_")

wtimages<-images[grep(x=images, "wt")]
leanimages<-images[grep(x=images, "lean_nnat")]                 
hugennatimages<-images[grep(x=images, "huge_nnat_")]
hugennatd9images<-images[grep(x=images, "huge_nnatd9")]

#wt images
wtimage1<-image_read(paste(p, wtimages[1], sep=""))
wtimage2<-c(wtimage1, image_read(paste(p, wtimages[2], sep="")))
wtimage3<-c(wtimage2, image_read(paste(p, wtimages[3], sep="")))
wtimage4<-c(wtimage3, image_read(paste(p, wtimages[4], sep="")))
wtimage5<-c(wtimage4, image_read(paste(p, wtimages[5], sep="")))
wtimage6<-c(wtimage5, image_read(paste(p, wtimages[6], sep="")))
finalwtimage<-image_append(wtimage6)
newimg<-image_border(image_background(finalwtimage, "white"), "#FFFFFF", "70x100")
finalwtimage<-image_annotate(newimg, "WT", size=72, gravity="North", color=colors[1], font="sans")


#lean images
leanimage1<-image_read(paste(p, leanimages[1], sep=""))
leanimage2<-c(leanimage1, image_read(paste(p, leanimages[2], sep="")))
leanimage3<-c(leanimage2, image_read(paste(p, leanimages[3], sep="")))
leanimage4<-c(leanimage3, image_read(paste(p, leanimages[4], sep="")))
finalleanimage<-image_append(leanimage4)
newimg<-image_border(image_background(finalleanimage, "white"), "#FFFFFF", "70x100")
finalleanimage<-image_annotate(newimg, "Lean", size=72, gravity="North", color=colors[2], font="sans")

#huge nnat images
hugennatimage1<-image_read(paste(p, hugennatimages[1], sep=""))
hugennatimage2<-c(hugennatimage1, image_read(paste(p, hugennatimages[2], sep="")))
hugennatimage3<-c(hugennatimage2, image_read(paste(p, hugennatimages[3], sep="")))
hugennatimage4<-c(hugennatimage3, image_read(paste(p, hugennatimages[4], sep="")))
finalhugennatimage<-image_append(hugennatimage4)
newimg<-image_border(image_background(finalhugennatimage, "white"), "#FFFFFF", "70x100")
finalhugennatimage<-image_annotate(newimg, "Huge Nnat", size=72, gravity="North", color=colors[3], font="sans")

#huge nnatd9 images
#hugennatd9image1<-image_read(paste(p, hugennatd9images[1], sep=""))
#hugennatd9image2<-c(hugennatd9image1, image_read(paste(p, hugennatd9images[2], sep="")))
#hugennatd9image3<-c(hugennatd9image2, image_read(paste(p, hugennatd9images[3], sep="")))
#hugennatd9image4<-c(hugennatd9image3, image_read(paste(p, hugennatd9images[4], sep="")))
#finalhugennatd9image<-image_append(hugennatd9image4)
#newimg<-image_border(image_background(finalhugennatd9image, "white"), "#FFFFFF", "70x100")
#finalhugennatd9image<-image_annotate(newimg, "Huge NnatD9", size=72, gravity="North", color=colors[3], font="sans")

promoterimage<-image_append(c(finalwtimage, finalleanimage, finalhugennatimage))#, finalhugennatd9image))
promoterimage<-image_border(image_background(promoterimage, "white"), "#FFFFFF", "70x150")
promoterimage<-image_annotate(promoterimage, paste(rt," Regions",sep=""), size=100, gravity="North", color=colors[1], font="sans")
image_write(promoterimage, path=paste(p,id,"_",rt,"_",hm,".png", sep=""), depth=16)


#enhancers
rm(list=ls())
library(magick)
colors=c("#000000", "#8299FF", "#0531FF")# wt, lean, giant colors
rt="Enhancer"
hm="H3K27ac"
dt="seqdepthnorm"
machine="local"
if (machine=="local"){
   p=paste("/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",hm,"_",dt,"/plots/", sep="")
}else{
   p=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",hm,"_",dt,"/plots/", sep="")
}
images<-list.files(path=p)

#filtering

images<-images[grep(x=images,"center")]
id=paste(strsplit(images, split="_")[[1]][3:4], collapse ="_")

wtimages<-images[grep(x=images, "wt")]
leanimages<-images[grep(x=images, "lean_nnat")]
hugennatimages<-images[grep(x=images, "huge_nnat_")]
hugennatd9images<-images[grep(x=images, "huge_nnatd9")]

#wt images
wtimage1<-image_read(paste(p, wtimages[1], sep=""))
wtimage2<-c(wtimage1, image_read(paste(p, wtimages[2], sep="")))
wtimage3<-c(wtimage2, image_read(paste(p, wtimages[3], sep="")))
wtimage4<-c(wtimage3, image_read(paste(p, wtimages[4], sep="")))
wtimage5<-c(wtimage4, image_read(paste(p, wtimages[5], sep="")))
wtimage6<-c(wtimage5, image_read(paste(p, wtimages[6], sep="")))
finalwtimage<-image_append(wtimage6)
newimg<-image_border(image_background(finalwtimage, "white"), "#FFFFFF", "70x100")
finalwtimage<-image_annotate(newimg, "WT", size=72, gravity="North", color=colors[1], font="sans")


#lean images
leanimage1<-image_read(paste(p, leanimages[1], sep=""))
leanimage2<-c(leanimage1, image_read(paste(p, leanimages[2], sep="")))
leanimage3<-c(leanimage2, image_read(paste(p, leanimages[3], sep="")))
leanimage4<-c(leanimage3, image_read(paste(p, leanimages[4], sep="")))
finalleanimage<-image_append(leanimage4)
newimg<-image_border(image_background(finalleanimage, "white"), "#FFFFFF", "70x100")
finalleanimage<-image_annotate(newimg, "Lean", size=72, gravity="North", color=colors[2], font="sans")

#huge nnat images
hugennatimage1<-image_read(paste(p, hugennatimages[1], sep=""))
hugennatimage2<-c(hugennatimage1, image_read(paste(p, hugennatimages[2], sep="")))
hugennatimage3<-c(hugennatimage2, image_read(paste(p, hugennatimages[3], sep="")))
hugennatimage4<-c(hugennatimage3, image_read(paste(p, hugennatimages[4], sep="")))
finalhugennatimage<-image_append(hugennatimage4)
newimg<-image_border(image_background(finalhugennatimage, "white"), "#FFFFFF", "70x100")
finalhugennatimage<-image_annotate(newimg, "Huge Nnat", size=72, gravity="North", color=colors[3], font="sans")

#huge nnatd9 images
#hugennatd9image1<-image_read(paste(p, hugennatd9images[1], sep=""))
#hugennatd9image2<-c(hugennatd9image1, image_read(paste(p, hugennatd9images[2], sep="")))
#hugennatd9image3<-c(hugennatd9image2, image_read(paste(p, hugennatd9images[3], sep="")))
#hugennatd9image4<-c(hugennatd9image3, image_read(paste(p, hugennatd9images[4], sep="")))
#finalhugennatd9image<-image_append(hugennatd9image4)
#newimg<-image_border(image_background(finalhugennatd9image, "white"), "#FFFFFF", "70x100")
#finalhugennatd9image<-image_annotate(newimg, "Huge NnatD9", size=72, gravity="North", color=colors[3], font="sans")

promoterimage<-image_append(c(finalwtimage, finalleanimage, finalhugennatimage))#, finalhugennatd9image))
promoterimage<-image_border(image_background(promoterimage, "white"), "#FFFFFF", "70x150")
promoterimage<-image_annotate(promoterimage, paste(rt," Regions",sep=""), size=100, gravity="North", color=colors[1], font="sans")
image_write(promoterimage, path=paste(p,id,"_",rt,"_",hm,".png", sep=""), depth=16)
            


