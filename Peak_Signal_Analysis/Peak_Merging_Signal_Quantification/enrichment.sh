#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=20
#PBS -N AUCs_2D
SECONDS=0

cd /secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/H3K27ac_seqdepthnorm/scripts
sh e${PBS_ARRAYID}.sh
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) %                  60))min          $(($SECONDS % 60))sec"
echo $ELAPSED
