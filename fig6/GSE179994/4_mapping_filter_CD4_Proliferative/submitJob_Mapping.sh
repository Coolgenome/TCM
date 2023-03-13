#BSUB -J submitJob_Mapping
#BSUB -q short
#BSUB -W 3:00
#BSUB -n 1
#BSUB -M 10
#BSUB -R rusage[mem=10]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE179994/4_mapping_filter_CD4_Proliferative/submitJob_Mapping.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE179994/4_mapping_filter_CD4_Proliferative/submitJob_Mapping.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE179994/4_mapping_filter_CD4_Proliferative/
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE179994/4_mapping_filter_CD4_Proliferative/submitJob_Mapping.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE179994/4_mapping_filter_CD4_Proliferative/submitJob_Mapping.e.txt
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
PIPELINE_NAME=GSE179994__4_mapping_filter_CD4_Proliferative
PIPELINE_PATH_NAME=GSE179994/4_mapping_filter_CD4_Proliferative
PROJECT_NAME=$(basename ${PROJECT_FOLDER})

OutDir=$RESULT_FOLDER/$PIPELINE_PATH_NAME
if [ ! -d $OutDir ]; then
    mkdir -p $OutDir
fi


runR="Rscript --no-save "

OUT_ROOT=$OutDir

qName=medium
wTime=24:00
cn=1
mem=100

dataPath=/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE179994/2_extractTcell_proliferative/CD4ProlifSeuratObj_2022-10-31.rds
dataFile=$(basename $dataPath)
extension="${dataFile##*.}"

referenceDataPath="/rsrch3/scratch/genomic_med/ychu2/projects/p1review/figureCode/result/0_write_sample_info/Proliferative_2022-10-20.rds"
JOBFOLDER=${OUT_ROOT}/proliferative
if [ ! -d $JOBFOLDER ]; then
    mkdir -p $JOBFOLDER
fi
JOBNAME=Mapping_proliferative
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
    /bin/bash -c "module load R/4.0.3; Rscript /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE179994/4_mapping_filter_CD4_Proliferative/Mapping.R -r ${referenceDataPath} -q ${dataPath} -o ${JOBFOLDER}"

