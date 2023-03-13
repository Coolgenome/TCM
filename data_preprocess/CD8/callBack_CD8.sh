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

rootFolder=$(dirname ${dataPath})


# ${runR} ${srcD}/iml_check_bubble-plot.R -d $dataPath -c "CD8_c0_Tcm CD8_c1_Tex CD8_c2_Teff CD8_c3_Naive CD8_c4_Stressed CD8_c5_ISG CD8_c6_Tcm_DKK3 CD8_c7_Pre-Tex CD8_c8_Teff_CX3CR1 CD8_c9_KLRC4 CD8_c10_Teff_CD244 CD8_c11_Teff_SEMA4A CD8_c12_Trm CD8_c13_Naive_TCF7 CD8_c14_gdT"

${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m  /rsrch3/home/genomic_med/ychu2/share/database/Pan-T/Reference/PMID34290406.txt

# tempMarkerFolder=/rsrch3/home/genomic_med/ychu2/configs/scSeqs/database/Markers/CD8/Genelistbyfunction
# for markerFile in $(ls $tempMarkerFolder); do
#     JOBNAME=job_${markerFile%%.*}
#     tempFolder=${rootFolder}
#     JOBFOLDER=${tempFolder}
#     if [ ! -d $tempFolder ]; then
#         mkdir -p $tempFolder
#     fi
#     if [ -f ${JOBFOLDER}/${JOBNAME}.o.txt ] || [ -f ${JOBFOLDER}/${JOBNAME}.e.txt ]; then
#         rm ${JOBFOLDER}/${JOBNAME}.*.txt -f
#     fi
#     bsub \
#             -J ${JOBNAME} \
#             -o ${JOBFOLDER}/${JOBNAME}.o.txt \
#             -e ${JOBFOLDER}/${JOBNAME}.e.txt \
#             -cwd ${JOBFOLDER} \
#             -q short \
#             -W 3:00 \
#             -n 1 \
#             -M 100\
#             -R rusage[mem=100] \
#             -B \
#             -N \
#             -u ychu2@mdanderson.org \
#             /bin/bash -c "module load python/3.7.3-anaconda; module load R/3.6.0; ${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m ${tempMarkerFolder}/${markerFile}; ${runR} ${srcD}/stack-ViolinPlot.R -d $dataPath -m ${tempMarkerFolder}/${markerFile}; ${runR} ${srcD}/feature-plot.R -d $dataPath -o $(dirname ${dataPath})/featureplot/${markerFile%%.*} -c ${paramD}/feature-plot-origin.json -m ${tempMarkerFolder}/${markerFile}"
# done

# ${runR} ${srcD}/marker-classification.R -d $(dirname ${dataPath})/snn-single-markers.tsv -m $tempMarkerFolder

# ${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m ${databaseD}/Pan-T/CD8/CD8Markers_byfunction_fig1A.txt
# ${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m ${databaseD}/Pan-T/CD8/CD8Markers_fig1A.txt
# ${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m ${databaseD}/Pan-T/CD8/CD8Markers.txt
# ${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m ${databaseD}/Markers/CD8/ProliferativeSignatures.txt

# ${runR} ${srcD}/stack-ViolinPlot.R -d $dataPath -m ${databaseD}/Pan-T/CD8/CD8Markers.txt


# tempMarkerFolder=/rsrch3/home/genomic_med/ychu2/configs/scSeqs/database/Markers/CD8/Genelistbyfunction

# for markerFile in $(ls $tempMarkerFolder); do
#     ${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m ${tempMarkerFolder}/${markerFile}
#     ${runR} ${srcD}/stack-ViolinPlot.R -d $dataPath -m ${tempMarkerFolder}/${markerFile}
#     ${runR} ${srcD}/feature-plot.R -d $dataPath -o $(dirname ${dataPath})/featureplot/${markerFile%%.*} -c ${paramD}/feature-plot-origin.json -m ${tempMarkerFolder}/${markerFile}
# done
