import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

hm="H3K27ac"
dt="seqdepthnorm"
df=pd.read_csv('/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/'+hm+"_"+dt+"/matrixfiles/enrichmentmat.tab", delimiter="\t")

#WT plot
dfnew=df.loc[:,['499_lean_wt', '551_lean_wt', '583_lean_wt', '674_lean_wt', '683_lean_wt','687_lean_wt']]

x=dfnew.mean(axis=1)
y=dfnew.std(axis=1)/x

x=x[x!=0]
y=y[~np.isnan(y)]
x=np.log10(x)
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
fig, ax = plt.subplots(figsize=(4,4))
#plotting
ax.scatter(np.log10(x), y, c=z, s=1, edgecolor='')
ax.set_xlabel('$log_{10}[mean(WT)]$')
ax.set_ylabel('$CV$')

#huge plot
dfnew=df.loc[:,['550_huge_nnat', '558_huge_nnat', '592_huge_nnat', '680_huge_nnat']]

x=dfnew.mean(axis=1)
y=dfnew.std(axis=1)/x

x=x[x!=0]
y=y[~np.isnan(y)]
x=np.log10(x)
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
fig, ax = plt.subplots(figsize=(4,4))
#plotting
ax.scatter(np.log10(x), y, c=z, s=1, edgecolor='')
ax.set_xlabel('$log_{10}[mean(Huge)]$')
ax.set_ylabel('$CV$')

#lean plot
dfnew=df.loc[:,['500_lean_nnat', '557_lean_nnat', '594_lean_nnat', '686_lean_nnat']]

x=dfnew.mean(axis=1)
y=dfnew.std(axis=1)/x

x=x[x!=0]
y=y[~np.isnan(y)]
x=np.log10(x)
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
fig, ax = plt.subplots(figsize=(4,4))
#plotting
ax.scatter(np.log10(x), y, c=z, s=1, edgecolor='')
ax.set_xlabel('$log_{10}[mean(Lean)]$')
ax.set_ylabel('$CV$')
