#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=28
#PBS -N average_plot_para

SECONDS=0

cd /secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/H3K27ac_seqdepthnorm

R --file="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/H3K27ac_seqdepthnorm/commonbed.R"
R --file="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/H3K27ac_seqdepthnorm/averageplot.R"

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) %                  60))min          $(($SECONDS % 60))sec"
echo $ELAPSED
