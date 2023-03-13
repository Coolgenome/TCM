#BSUB -J harmony
#BSUB -q e80medium
#BSUB -W 23:50
#BSUB -n 1
#BSUB -M 550
#BSUB -R rusage[mem=550]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT2_harmony/1_run_harmony/harmony.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT2_harmony/1_run_harmony/harmony.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT2_harmony/1_run_harmony/
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT2_harmony/1_run_harmony/harmony.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT2_harmony/1_run_harmony/harmony.e.txt
module load python/3.7.3-anaconda
module load R/4.0.3
#____----____----____

PROJECT_FOLDER=/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7
DATA_FOLDER=${PROJECT_FOLDER}/data
RESULT_FOLDER=${PROJECT_FOLDER}/result
CODE_FOLDER=${PROJECT_FOLDER}/code
PIPELINE_FOLDER=${CODE_FOLDER}/pipeline
SRC_FOLDER=${CODE_FOLDER}/src
KNOWLEDGE_FOLDER=${PROJECT_FOLDER}/knowledge
PIPELINE_NAME=GSE169246__subT2_harmony__1_run_harmony
PIPELINE_PATH_NAME=GSE169246/subT2_harmony/1_run_harmony
PROJECT_NAME=$(basename ${PROJECT_FOLDER})

OutDir=$RESULT_FOLDER/$PIPELINE_PATH_NAME
if [ ! -d $OutDir ]; then
    mkdir -p $OutDir
fi


Rscript --no-save /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT2_harmony/1_run_harmony/harmony.R \
        -d /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/subT1_extract_from_b2/1_extract/outs/subObj_2022-09-01.rds \
        -o $OutDir/harmony_output.rds
