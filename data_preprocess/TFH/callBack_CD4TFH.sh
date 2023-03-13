##!/usr/bin/env bash

projectPath=/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject
analysisPath=${projectPath}/analysis
mainscriptsFolder=${analysisPath}/scripts
pipelinesFolder=${mainscriptsFolder}/pipelines
srcD=${analysisPath}/scripts/src
paramD=${analysisPath}/scripts/params
databaseD='/rsrch3/home/genomic_med/ychu2/share/database'

runR="Rscript --no-save "

dataPath=${1}

# ${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m ${databaseD}/THelper.txt
# ${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m ${databaseD}/TCD4.txt
# ${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m ${databaseD}/Pan-T/CD4/TFH/CD4TFHMarkers.txt
# ${runR} ${srcD}/stack-ViolinPlot.R -d $dataPath -m ${databaseD}/Pan-T/CD4/TFH/CD4TFHMarkers.txt

# JOBNAME=job_monocle3
# tempFolder=$(dirname $dataPath)
# JOBFOLDER=${tempFolder}
# if [ ! -d $tempFolder ]; then
#     mkdir -p $tempFolder
# fi
# if [ -f ${JOBFOLDER}/${JOBNAME}.o.txt ] || [ -f ${JOBFOLDER}/${JOBNAME}.e.txt ]; then
#     rm ${JOBFOLDER}/${JOBNAME}.*.txt -f
# fi
# bsub \
#     -J ${JOBNAME} \
#     -o ${JOBFOLDER}/${JOBNAME}.o.txt \
#     -e ${JOBFOLDER}/${JOBNAME}.e.txt \
#     -cwd ${JOBFOLDER} \
#     -q long \
#     -W 120:00 \
#     -n 1 \
#     -M 500\
#     -R rusage[mem=500] \
#     -B \
#     -N \
#     -u ychu2@mdanderson.org \
#     /bin/bash -c "module load python/3.7.3-anaconda; module load R/3.6.0; ${runR} ${srcD}/trajectory_monocle3_from_seurat.r -d ${dataPath} -c 3"
