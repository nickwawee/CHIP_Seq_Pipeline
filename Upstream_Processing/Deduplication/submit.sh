#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=1
#PBS -N dedup

SECONDS=0

cd /home/nick.wawee/nickS/Leslie.Yang/adipocyte_data/deduplication2/scripts/
sh dd${PBS_ARRAYID}.sh
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) %                  60))min          $(($SECONDS % 60))sec"
echo $ELAPSED
