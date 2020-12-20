#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=4
#PBS -N ChIP-seq

SECONDS=0

cd /secondary/projects/mnp/nick/testing/relacs/relacs/

python demultiplex_relacs.py --umiLength 4 sampleTable.tsv output

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) %                  60))min          $(($SECONDS % 60))sec"
echo $ELAPSED

