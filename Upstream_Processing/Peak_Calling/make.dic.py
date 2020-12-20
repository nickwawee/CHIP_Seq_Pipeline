import os,sys
import argparse

ap=argparse.ArguementParser()
ap.add_argument("-bd", "--bamdir", required=True, help="directory with bam files that have the name in the following format, animalnumber_phenotype_genotype_sex_antibody")
ap.add_argument("-od", "--outputdir", required=True, help="directory that the yaml file goes to")
args = vars(ap.parse_args())

x=[f for f in os.listdir(str(args['bamdir'])) if f.endswith("bam")]
x.sort()

if str(args['outputdir']).endswith("/"):
   outdir=str(args['outputdir'])
else:
   outdir=str(args['outputdir'])+"/" 

d={}
for f in x:
	f=f[:-19]
	prefix=f.split("_")[0]
	if prefix in d:
		L=d[prefix]
		L.append(f)
		d[prefix]=L
	else:
		d[prefix]=[f]	
h=open(outdir+"ChIPdict.yaml", "w")
h.write("chip_dict:\n")

for k in d.keys():
	L=d[k]
	Input= [s for s in L if "_Input" in s]
	if len(Input)>1:
		sys.exit()
	for s in L:
		if "_Input" not in s:
			h.write("   "+s+":\n")
			h.write("      "+"control: "+Input[0]+"\n")
			h.write("      "+"broad: False\n")

h.close()
