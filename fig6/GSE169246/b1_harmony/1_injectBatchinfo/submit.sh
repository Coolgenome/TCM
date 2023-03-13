#BSUB -J submit
#BSUB -q e80medium
#BSUB -W 23:00
#BSUB -n 1
#BSUB -M 550
#BSUB -R rusage[mem=550]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/b1_harmony/1_injectBatchinfo/submit.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/b1_harmony/1_injectBatchinfo/submit.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/b1_harmony/1_injectBatchinfo/
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/b1_harmony/1_injectBatchinfo/submit.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/b1_harmony/1_injectBatchinfo/submit.e.txt
module load python/3.7.3-anaconda
module load R/4.0.3
#____----____----____


Rscript --no-save /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/b1_harmony/1_injectBatchinfo/inject.R
