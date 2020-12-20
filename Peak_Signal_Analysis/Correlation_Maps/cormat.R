#This script will create a correlation matrix between samples based on their AUC's
#http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
rm(list=ls())
library(reshape2)
library(colorRamps)
library(ggplot2)
load("/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/matrixfiles/enrichmentdf.RData")
hm="H3K27ac"
load("/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/mousedata.Rdata")
mousedata=features
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd, method='ward.D2')
  cormat <-cormat[hc$order, hc$order]
}

#tidying

enrichmentdf<-enrichmentdf[,-ncol(enrichmentdf)]
nums<-as.numeric(unlist(lapply(1:ncol(enrichmentdf), function(i){
  strsplit(colnames(enrichmentdf), split="_")[[i]][1]
})))

order="pheno"
if (order=="pheno"){
  enrichmentdf<-enrichmentdf[, match(mousedata$mousenumber, nums)]
  #making matrix
  cormat<-round(cor(enrichmentdf),2)
  # Reorder the correlation matrix
  #cormat <- reorder_cormat(cormat)#use this when sorting via ward.D2
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradientn(colours = colorRamps::matlab.like(100), limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1), axis.text.y = element_text(size=10))+
    coord_fixed()+
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))+
    ggtitle(paste(hm,' Correlation Map, Peak AUCs', sep=""))

# Print the heatmap
  print(ggheatmap)
  ggsave(filename=paste("/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/plots/correlationmatrix_",order,"_",hm,".png", sep=""), units = "in", width = 7, height = 6, dpi = 600)
}

order="ward"
if (order=="ward"){
  enrichmentdf<-enrichmentdf[, match(mousedata$mousenumber, nums)]
  #making matrix
  cormat<-round(cor(enrichmentdf),2)
  # Reorder the correlation matrix
  cormat <- reorder_cormat(cormat)#use this when sorting via ward.D2
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradientn(colours = colorRamps::matlab.like(100), limit = c(0,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 10, hjust = 1), axis.text.y = element_text(size=10))+
    coord_fixed()+
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))+
    ggtitle(paste(hm,' Correlation Map, Peak AUCs', sep=""))
  
  # Print the heatmap
  print(ggheatmap)
  ggsave(filename=paste("/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq/comp_mat/bedregions_mat/plots/correlationmatrix_",order,"_",hm,".png", sep=""), units = "in", width = 7, height = 6, dpi = 600)
}