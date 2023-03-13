#BSUB -J merge
#BSUB -q highmem
#BSUB -W 23:00
#BSUB -n 1
#BSUB -M 100
#BSUB -R rusage[mem=100]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE144649/1_merge/merge.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE144649/1_merge/merge.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE144649/1_merge/
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE144649/1_merge/merge.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE144649/1_merge/merge.e.txt
module load python/3.7.3-anaconda
module load R/3.6.0
#____----____----____

OutFolder=/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE144649/1_merge
if [ ! -d $OutFolder ]; then
    mkdir -p $OutFolder
fi

Rscript --no-save  /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE144649/1_merge/merge.R -d /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/data/GSE144469/raw -o $OutFolder
