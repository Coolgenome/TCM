##!/usr/bin/env bash

module load python/3.7.3-anaconda
module load R/4.0.3

mainscriptFolder=${1}
inData=${2}
currentFolder=${3}
dist=${4}
reduction=${5}
npc=${6}
parentJobName=${7}
ClusterArray=${8}

IFS=';' read -a ClusterArray <<< "$ClusterArray"

nneighbors=${9}
toRunUMAP=${10}
toRunClustering=${11}
toRunCommonAnalysis=${12}
toRunCallBack=${13}
callBackPath=${14}

echo "RunUMAPJobs parameters:
mainscriptFolder=${mainscriptFolder}
inData=${inData}
currentFolder=${currentFolder}
dist=${dist}
reduction=${reduction}
npc=${npc}
parentJobName=${parentJobName}
ClusterArray[@]=${ClusterArray[@]}
nneighbors=${nneighbors}
toRunUMAP=${toRunUMAP}
toRunClustering=${toRunClustering}
toRunCommonAnalysis=${toRunCommonAnalysis}
toRunCallBack=${toRunCallBack}
callBackPath=${callBackPath}
"

srcD=${HOME}/configs/public/pipeline/UMAP_CLUSTER_JOBS_EMBEDED
runR="Rscript --no-save "

if [ $toRunUMAP = "YES" ]; then
    ${runR} ${srcD}/RunUMAP.R -d ${inData} -o ${currentFolder}/umap.rds -r ${reduction} -n ${npc} -i ${dist} -e ${nneighbors}
fi

inData=${currentFolder}/umap.rds
mmForInData=`du --apparent-size --block-size=1000000000 ${inData}  | awk '{print $1}'` # GB
if [ "$mmForInData" -gt "1200" ]; then
    mmForInData=350
fi

hoursTime=$(( ${mmForInData} + 24 ))
wTime="${hoursTime}:00"

qName="e80short"
if [ "${hoursTime}" -lt "2" ]; then
    qName="e80short"
else if [ "${hoursTime}" -lt "24" ]; then
         qName="e80medium"
     else if [ "${hoursTime}" -lt "120" ]; then
              qName="highmem"
          else
              qName="vhighmem"
          fi
     fi
fi

if [ "$toRunClustering" == "NO" ]; then
    mmForInData=320
    qName="highmem"
    wTime="72:00"
fi

rootFolder=$currentFolder
for res in "${ClusterArray[@]}"; do
    JOBNAME=${parentJobName}_CLUSTER_res_${res}
    tempFolder=${rootFolder}/${JOBNAME}
    if [ ! -d $tempFolder ]; then
        mkdir -p $tempFolder
    fi
    JOBFOLDER=${tempFolder}
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
        -n 1 \
        -M $(( ${mmForInData} + 200 ))\
        -R rusage[mem=$(( ${mmForInData} + 200 ))] \
        -B \
        -N \
        -u ychu2@mdanderson.org \
        /bin/bash -c "${HOME}/configs/public/pipeline/UMAP_CLUSTER_JOBS_EMBEDED/FindClusterJobs.sh ${mainscriptFolder} ${inData} ${tempFolder} ${res} ${reduction} ${npc} ${JOBNAME} ${toRunClustering} ${toRunCommonAnalysis} ${toRunCallBack} ${callBackPath}"
done
