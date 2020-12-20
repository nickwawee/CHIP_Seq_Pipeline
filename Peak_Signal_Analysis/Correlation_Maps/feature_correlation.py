#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 10:51:31 2019

@author: nick.wawee

This script will load correlation columns from a given directory and perform kmeans clustering and then hierarchialclusterings on centriods if desired.
Its input is a tab separated matrix with features as rows and samples as columns. Its output is a heatmap that is hierarchial clustered  by the Ward agglomerative clustering method. It takes about 8 hours for a 24000x18 matrix without using the kmeans clustering. An improvement that could be made is adding the silhouette method.
"""

import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


hm="H3K27ac"
dt="seqdepthnorm"

df=pd.read_csv("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/matrixfiles/enrichmentmat_filtvar.tab", sep="\t")


#df2=df.sample(50000) #use when testing
df3=df.loc[(df.abs().sum(axis=1) != 0),:]#removing all 0 regions

corrmat=df3.transpose().corr(method='spearman')

#1st heatmap
heatmap=sns.clustermap(corrmat, metric="euclidean", method="ward")
heatmap.savefig("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/plots/correlation.png",dpi=600)
'''
#use if hierarchial does not work, finding elbow
numclust=20
ssd=[]
for i in range(1,numclust):
    kmeans=KMeans(n_clusters=i, random_state=123)
    ssd.append(kmeans.fit(corrmat).inertia_)

plt.scatter(range(1,numclust), ssd)
plt.xlabel('k')
plt.ylabel('Total Within Variation Sum of Squares')
plt.xticks(range(1,numclust,2))
plt.savefig("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/plots/correlation_elbow.png",dpi=600)

#fitting kmeans
numclust=6
kmeans=KMeans(n_clusters=numclust, random_state=123)
res=kmeans.fit(corrmat)
centriods=res.cluster_centers_
ssd=res.inertia_

#hierarchial clustering on centriods
cluster = AgglomerativeClustering(n_clusters=numclust, affinity='euclidean', linkage='ward')
clusters=cluster.fit_predict(centriods.transpose())
clustersdf=pd.DataFrame(clusters)
clustersdf.index=corrmat.index

#sorting rows
alldf=pd.concat([clustersdf,corrmat], axis=1)
alldf=alldf.sort_values(by=[0])
alldf=alldf.drop([0], axis=1)

#sorting columns
alldf=pd.concat([clustersdf.transpose(),alldf])
alldf=alldf.sort_values(by=[0], axis=1)
alldf=alldf.drop([0], axis=0)

#plotting heatmap
heatmap=sns.clustermap(alldf.abs(), row_cluster=False, col_cluster=False)
#clustersdf=clustersdf.sort_values(by=[0])
heatmap.savefig("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/"+hm+"_"+dt+"/plots/correlation.png",dpi=600)
'''
