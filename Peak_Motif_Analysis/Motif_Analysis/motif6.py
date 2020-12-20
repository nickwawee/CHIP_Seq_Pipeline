import os, sys
from Bio import SeqIO
from subprocess import *
import argparse

#getting the sequences
def get_sequence(bedReg,genomeFas,outDir,regType,bp):

    seqs=[s for s in SeqIO.parse(open(genomeFas),"fasta")]
    lines=open(bedReg).readlines()
    h=open(outDir+bedReg.split("/")[-1].split(".bed")[0] +".fas","w")
    count=0
    for l in lines:
        if l.startswith("chr"):
            cr=l.split("\t")[0]
            st=int(l.split("\t")[1])
            en=int(l.split("\t")[2])
            for seq in seqs:
                if seq.id==cr:
                    s=str(seq.seq)
                    if regType=="center":
                        cen=int((st+en)/2)
                        if (cen> bp+1) and cen < (len(seq)-(bp+1)):
                            h.write(">seq_"+'{:06}'.format(count)+"\n")
                            h.write(s[(cen-bp):(cen+bp)]+"\n")
                            count+=1
                    elif regType=="whole":
                        h.write(">seq_"+'{:06}'.format(count)+"\n")
                        h.write(s[st:en]+"\n")
                        count+=1
    h.close()



#Runing MEME - AME
def runMeme(fasReg,nc,nm,control,wd):#fasReg=region sequence files, number of cores, number of motifs, sequences to be compared against, working directory for meme

    h=open(wd+"meme.run.sh","w")
    h.write("source /primary/projects/mnp/tools/anaconda3/etc/profile.d/conda.sh\n")
    h.write("conda activate meme_env\n")
    h.write("cd "+wd+"\n")
    h.write("meme "+fasReg+" -dna "+ "-p "+str(nc)+" -mod anr "+"-nmotifs "+str(nm)+" -objfun classic "+"-o "+wd+"memeoutput"+"\n" )#meme command, meme giant.fas.masked -dna -p 28 -mod anr -nmotifs 20 -objfun de -neg wt.fas.masked -o memeoutput_giant_diffenrich
    h.write("ame --o "+wd+"ameoutput"+" --control "+control+" "+fasReg+" "+ wd+"memeoutput/meme.txt")#ame command, ame --o ameoutput --control wt.fas.masked giant.fas.masked memeoutput_giant/meme.txt
    h.close()
    cmd =["/bin/bash", wd+"meme.run.sh"]

    proc=Popen(cmd,stdout=PIPE,stderr=PIPE)
    (output, error)=proc.communicate()
    return_code = proc.wait()
    if return_code != 0:
                sys.stderr.write(error)
                sys.exit(1)

#Running HOMER
def runHomer(Reg,control,genomeFas,outDir,preparseDir,filetype,nc,nm):#input bed/fasta regions, genome fasta file, output directory, preparsed directory, bed or fasta file type, num of cores, num of motifs
    h=open(outDir+"homer.run.sh", "w")
    h.write("source /primary/projects/mnp/tools/anaconda3/etc/profile.d/conda.sh\n")
    h.write("conda activate homer_env\n")
    if filetype=="bed":
        h.write("findMotifsGenome.pl "+Reg+" "+genomeFas+" "+outDir+" -size "+" given "+ " -preparsedDir "+preparseDir+ " -p "+str(nc))
        h.close()
        cmd =["/bin/bash", outDir+"homer.run.sh"]
        proc=Popen(cmd,stdout=PIPE,stderr=PIPE)
        (output, error)=proc.communicate()
        return_code = proc.wait()
        if return_code != 0:
                 sys.stderr.write(error)
                 sys.exit(1)
    if filetype=="fasta":
        h.write("findMotifs.pl "+Reg+" fasta "+ outDir + " -fasta "+control+" -S "+str(nm)+" -p "+str(nc)+" -nomask ")
        h.close()
        cmd =["/bin/bash", outDir+"homer.run.sh"]
        proc=Popen(cmd,stdout=PIPE,stderr=PIPE)
        (output, error)=proc.communicate()
        return_code = proc.wait()
        if return_code != 0:
                 sys.stderr.write(error)
                 sys.exit(1)

#Running RepeatMasker, returns Ns for repeats
def runRepeatMasker(fasReg,outDir,nc,species):#input Fasta ,output directory, number of cores, species 

    cmd=["RepeatMasker", fasReg, "-species", species, "-pa", str(nc),"-dir", outDir]
    proc=Popen(cmd,stdout=PIPE,stderr=PIPE)
    (output, error)=proc.communicate()
    return_code = proc.wait()
    if return_code != 0:
                 sys.stderr.write(error)
                 sys.exit(1)


# Construct the argument parser
ap = argparse.ArgumentParser()

# Add the arguments to the parser
ap.add_argument("-p", "--processors", required=True,
   help="number of processors")
ap.add_argument("-m", "--motifs", required=True,
   help="number of motifs")
ap.add_argument("-b", "--bedfiles", required=True,
   help="bedfiles directory, may only contain two bed files of interest, the primary sequences and control. Files must be named primary.bed and control.bed. Must have / on end.")
ap.add_argument("-g", "--genome", required=True,
   help="genome.fa file")
ap.add_argument("-o", "--output", required=True,
   help="output directory, must have / on end.")
ap.add_argument("-s", "--size", required=True,
   help="size of region in basepairs centered on the peak center. Ex.) size of 1000 would return sequences of 1000 bp in length which are centered on the peak center. Enter 0 for whole sequences to be analyzed.")
args = vars(ap.parse_args())

#define arguments
nc=int(args['processors'])
nm=int(args['motifs'])
bedReg_=str(args['bedfiles'])
genomeFas_=str(args['genome'])
outDir_=str(args['output'])
size=int(args['size'])
#makes directories for homer and meme
cmd =["mkdir", outDir_+"homer"]
proc=Popen(cmd,stdout=PIPE,stderr=PIPE)
(output, error)=proc.communicate()
return_code = proc.wait()
if return_code != 0:
    sys.stderr.write(error)
    sys.exit(1)
cmd =["mkdir", outDir_+"meme"]
proc=Popen(cmd,stdout=PIPE,stderr=PIPE)
(output, error)=proc.communicate()
return_code = proc.wait()
if return_code != 0:
    sys.stderr.write(error)
    sys.exit(1)

#getting sequences and running repeatmasker
if size > 0:
    get_sequence(bedReg_+"primary.bed",genomeFas_,outDir_,"center",int(size/2))
    runRepeatMasker(outDir_+"primary.fas",outDir_,nc,"mouse")
    print('primary masking complete')
    get_sequence(bedReg_+"control.bed",genomeFas_,outDir_,"center",int(size/2))
    runRepeatMasker(outDir_+"control.fas",outDir_,nc,"mouse")
    print('control masking complete')
if size == 0:
    get_sequence(bedReg_+"primary.bed",genomeFas_,outDir_,"whole",size)
    runRepeatMasker(outDir_+"primary.fas",outDir_,nc,"mouse")
    print('primary masking complete')
    get_sequence(bedReg_+"control.bed",genomeFas_,outDir_,"whole",size)
    runRepeatMasker(outDir_+"control.fas",outDir_,nc,"mouse")
    print('control masking complete')
#running Meme and Homer
runMeme(outDir_+"primary.fas.masked", nc, nm, outDir_+"control.fas.masked", outDir_+"meme/")#(fasReg,nc,nm,control,wd)
print('Meme job complete')
runHomer(outDir_+"primary.fas.masked",outDir_+"control.fas.masked","NA",outDir_+"homer/",outDir_+"homer/preparsed/","fasta",nc,nm)
print('Homer job complete, end of job')
