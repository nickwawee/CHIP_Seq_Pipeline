#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=1
#PBS -N merge

SECONDS=0

cd /secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/relacs/combined2/scripts2

sh m${PBS_ARRAYID}.sh
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) %                  60))min          $(($SECONDS % 60))sec"
echo $ELAPSED
(base) [nick.wawee@submit WGB-seq]$

