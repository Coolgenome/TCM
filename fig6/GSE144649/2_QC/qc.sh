#BSUB -J qc
#BSUB -q short
#BSUB -W 2:30
#BSUB -n 1
#BSUB -M 100
#BSUB -R rusage[mem=100]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE144649/2_QC/qc.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE144649/2_QC/qc.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE144649/2_QC/
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE144649/2_QC/qc.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE144649/2_QC/qc.e.txt
module load python/3.7.3-anaconda
module load R/3.6.0
#____----____----____

OutDir=/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE144649/2_QC
if [ ! -d $OutDir ]; then
    mkdir -p $OutDir
fi

Rscript --no-save  /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/src/public/filter-QC.R \
        -d /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE144649/1_merge/merged.rds \
        -o $OutDir/qc.rds
