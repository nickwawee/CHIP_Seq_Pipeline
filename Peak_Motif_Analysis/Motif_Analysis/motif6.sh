#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=40
#PBS -N motif_6
#PBS -q longq
cd /secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/H3K27ac_PeakCompare_Motif/
SECONDS=0
python motif6.py -p 40 -m 100 -g /primary/projects/mnp/tools/data/repository/organisms/GRCm38.p6_GENCODE/genome_fasta/genome.fa -b bedfiles/giant_leanwt/ -o /secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/H3K27ac_PeakCompare_Motif/ -s 0
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) %                  60))min          $(($SECONDS % 60))sec"
echo $ELAPSED
