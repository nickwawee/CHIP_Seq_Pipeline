import os
import argparse

counter=0

ap = argparse.ArgumentParser()

# Add the arguments to the parser
ap.add_argument("-id", "--inputdir", required=True,
   help="input directory that contains bam files")
ap.add_argument("-od", "--outputdir", required=True,
   help="output directory for deduplicated bam files")
ap.add_argument("-sp", "--scriptdir", required=True,
   help="output directory for scripts to be submitted to PBS")
args = vars(ap.parse_args())

#define arguments
p1=str(args['inputdir'])
op=str(args['outputdir'])
sh_p=str(args['scriptdir'])

#retrieves all input file names
fq1=[r for r in os.listdir(p1) if r.endswith(".bam")]
for name in fq1:
    f=open(sh_p+"dd"+str(counter)+".sh","w")
    f.write("umi_tools dedup --paired -I "+p1+name +" -S "+op+name.split(".bam")[0]+"_dedup.bam")
    f.close()
    counter+=1
