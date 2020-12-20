
import os

#sample folder names, required to identify unmerged files
L1=["Sample_18L009599","Sample_18L009600","Sample_18L009601","Sample_18L009602","Sample_18L009603"]
L2=["Sample_18L010651","Sample_18L010652","Sample_18L010653","Sample_18L010654","Sample_18L010655"]
L3=["Sample_18L010651","Sample_18L010652","Sample_18L010653","Sample_18L010654","Sample_18L010655"] 
L4=["Sample_18L010651","Sample_18L010652","Sample_18L010653","Sample_18L010654","Sample_18L010655"]


counter=0
for j in range(len(L4)):
    p1="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/relacs/181031/output/"+L1[j]
    fq1=[r for r in os.listdir(p1) if r.endswith("fastq.gz")]
    
    p2="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/relacs/190116/output/"+L2[j]
    fq2=[r for r in os.listdir(p2) if r.endswith("fastq.gz")]
    
    p3="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/relacs/190123/output/"+L3[j]
    fq3=[r for r in os.listdir(p3) if r.endswith("fastq.gz")]
    
    p4="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/relacs/190125/output/"+L4[j]
    fq4=[r for r in os.listdir(p4) if r.endswith("fastq.gz")]
    
    op="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/relacs/combined2/"#merged files output directory 
    sh_p="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/relacs/combined2/scripts2/"#output script directory
    
    for name in fq1:
        ip1="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/relacs/181031/output/"+L1[j]+"/"
        ip2="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/relacs/190116/output/"+L2[j]+"/"
        ip3="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/relacs/190123/output/"+L3[j]+"/"
        ip4="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/relacs/190125/output/"+L4[j]+"/"
        
        f=open(sh_p+"m"+str(counter)+".sh","w")
        f.write("zcat  "+ip1+name + " "+ip2+name+ " "+ip3+name+ " "+ip4+name + "  | gzip  -n9  > " +op+name)
        f.close()
        counter+=1
