#BSUB -J submitJob_Mapping
#BSUB -q short
#BSUB -W 1:00
#BSUB -n 1
#BSUB -M 10
#BSUB -R rusage[mem=10]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT3_mapping_filter_split_by_marker/submitJob_Mapping.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT3_mapping_filter_split_by_marker/submitJob_Mapping.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT3_mapping_filter_split_by_marker/
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT3_mapping_filter_split_by_marker/submitJob_Mapping.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT3_mapping_filter_split_by_marker/submitJob_Mapping.e.txt
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
PIPELINE_NAME=GSE169246__subT3_mapping_filter_split_by_marker
PIPELINE_PATH_NAME=GSE169246/subT3_mapping_filter_split_by_marker
PROJECT_NAME=$(basename ${PROJECT_FOLDER})

OutDir=$RESULT_FOLDER/$PIPELINE_PATH_NAME
if [ ! -d $OutDir ]; then
    mkdir -p $OutDir
fi

runR="Rscript --no-save "

QUERYDATAFOLDER="/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/subT2_split_by_marker/outs"

CD8_ReferenceDataPath="/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/CD8_V6/nPC_50/UMAP_dist_0.1_nneighbor_50/p1CD8_V6_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds"
CD4_ReferenceDataPath="/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/CD4_V7/nPC_50/UMAP_dist_0.1_nneighbor_50/p1CD4_V7_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds"

PIPELINE_FOLDER="/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT3_mapping_filter_split_by_marker"

qName=e80long
wTime=100:00
cn=1
mem=300

for dataPath in ${QUERYDATAFOLDER}/*; do
    dataFile=$(basename $dataPath)
    extension="${dataFile##*.}"
    filename="${dataFile%.*}"
    referenceDataPath=""
    if [ "${filename}" = "CD8" ]; then
        referenceDataPath=${CD8_ReferenceDataPath}
    fi
    if [ "${filename}" = "CD4" ]; then
        referenceDataPath=${CD4_ReferenceDataPath}
    fi
    JOBFOLDER=${OutDir}/${filename}
    if [ ! -d $JOBFOLDER ]; then
        mkdir -p $JOBFOLDER
    fi

    JOBNAME=Mapping_filter_${filename}
    if [ -f ${JOBFOLDER}/${JOBNAME}.o.txt ] || [ -f ${JOBFOLDER}/${JOBNAME}.e.txt ]; then
        rm ${JOBFOLDER}/${JOBNAME}.*.txt -f
    fi

    bsub \
        -J ${JOBNAME} \
        -o ${JOBFOLDER}/${JOBNAME}.o.txt \
        -e ${JOBFOLDER}/${JOBNAME}.e.txt \
        -cwd ${JOBFOLDER} \
        -q ${qName} \
        -W ${wTime} \
        -n ${cn} \
        -M ${mem} \
        -R rusage[mem=${mem}] \
        -B \
        -N \
        -u ychu2@mdanderson.org \
        /bin/bash -c "module load R/4.0.3; Rscript ${PIPELINE_FOLDER}/Mapping.R -r ${referenceDataPath} -q ${dataPath} -o ${JOBFOLDER}"
done

