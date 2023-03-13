#BSUB -J pca
#BSUB -q short
#BSUB -W 2:50
#BSUB -n 1
#BSUB -M 100
#BSUB -R rusage[mem=100]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/3_pca/pca.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/3_pca/pca.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/3_pca/
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/3_pca/pca.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/3_pca/pca.e.txt
module load python/3.7.3-anaconda
module load R/3.6.0
#____----____----____


OutDir=/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/3_pca
if [ ! -d $OutDir ]; then
    mkdir -p $OutDir
fi

Rscript --no-save /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/3_pca/pca.R \
        -d /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/2_QC/qc.rds \
        -o $OutDir
