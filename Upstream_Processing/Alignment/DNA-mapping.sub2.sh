#PBS -l walltime=96:00:00
#PBS -l nodes=1:ppn=40
#PBS -N DNA-mapping2

cd /home/nick.wawee/nickS/Leslie.Yang/adipocyte_data/DNAmapping2/

SECONDS=0


DNA-mapping -i /home/nick.wawee/nickS/Leslie.Yang/adipocyte_data/DNAmapping/input  -o /home/nick.wawee/nickS/Leslie.Yang/adipocyte_data/DNAmapping2/output  --mapq 5 -j 40 mm10 --local --fastqc


ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) %                  60))min          $(($SECONDS % 60))sec"
 echo $ELAPSED
