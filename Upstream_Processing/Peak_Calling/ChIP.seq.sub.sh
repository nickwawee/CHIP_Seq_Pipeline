#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=28
#PBS -N ChIP-seq

cd /home/nick.wawee/nickS/Leslie.Yang/adipocyte_data/DNAmapping2/output

SECONDS=0

ChIP-seq -d /home/nick.wawee/nickS/Leslie.Yang/adipocyte_data/DNAmapping2/output mm10 /home/nick.wawee/nickS/Leslie.Yang/adipocyte_data/DNAmapping2/output/ChIPdict.yaml --local -j 28 --fromBam /secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/deduplication2/output/ --fromBamExt '.filtered_dedup.bam'


ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) %                  60))min          $(($SECONDS % 60))sec"
echo $ELAPSED
