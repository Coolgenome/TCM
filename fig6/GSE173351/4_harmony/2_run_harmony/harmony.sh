#BSUB -J harmony
#BSUB -q highmem
#BSUB -W 23:00
#BSUB -n 1
#BSUB -M 100
#BSUB -R rusage[mem=100]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/4_harmony/2_run_harmony/harmony.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/4_harmony/2_run_harmony/harmony.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/4_harmony/2_run_harmony/
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/4_harmony/2_run_harmony/harmony.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/4_harmony/2_run_harmony/harmony.e.txt
module load python/3.7.3-anaconda
module load R/3.6.0
#____----____----____

OutDir=/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/4_harmony
if [ ! -d $OutDir ]; then
    mkdir -p $OutDir
fi

Rscript --no-save /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/src/public/run-harmony.R \
        -d /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/4_harmony/harmony_input.rds \
        -o $OutDir/harmony_output.rds
