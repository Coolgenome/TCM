##!/usr/bin/env bash


srcD=/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/src/public

runR="Rscript --no-save "

dataPath=${1}

echo "begin"
${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m  /rsrch3/home/genomic_med/ychu2/share/database/Pan-T/Reference/PMID34290406.txt

# rootFolder=$(dirname ${dataPath})
# tempMarkerFolder=/rsrch3/home/genomic_med/ychu2/configs/scSeqs/database/Markers/CD4/Genelistbyfunction
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

# ${runR} ${srcD}/marker-classification.R -d $(dirname ${dataPath})/snn-single-markers.tsv -m /rsrch3/home/genomic_med/ychu2/configs/scSeqs/database/Markers/CD4/Genelistbyfunction

# ${runR} ${srcD}/tissue-composition-plot.R -d $dataPath
# ${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m ${databaseD}/Markers/CD4/markers/CD4_naive_clusters_comparison.txt

# ${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m ${databaseD}/TCD4.txt
# ${runR} ${srcD}/stack-ViolinPlot.R -d $dataPath -m ${databaseD}/Pan-T/CD4/CD4Markers.txt
# ${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m ${databaseD}/Pan-T/CD4/CD4Markers.txt
# ${runR} ${srcD}/stack-ViolinPlot.R -d $dataPath -m ${databaseD}/Pan-T/CD4/CD4Markers.txt
# ${runR} ${srcD}/monocleForTest.R -d $dataPath -n 8
