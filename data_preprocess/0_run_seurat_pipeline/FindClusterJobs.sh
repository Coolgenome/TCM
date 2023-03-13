##!/usr/bin/env bash

module load python/3.7.3-anaconda
module load R/4.0.3

mainscriptFolder=${1}
inData=${2}
currentFolder=${3}
res=${4}
reduction=${5}
npc=${6}
parentJobName=${7}
toRunClustering=${8}
toRunCommonAnalysis=${9}
toRunCallBack=${10}
callBackPath=${11}

echo "FindClusterJobs parameters:
mainscriptFolder=${mainscriptFolder}
inData=${inData}
currentFolder=${currentFolder}
res=${res}
reduction=${reduction}
npc=${npc}
parentJobName=${parentJobName}
toRunClustering=${toRunClustering}
toRunCommonAnalysis=${toRunCommonAnalysis}
toRunCallBack=${toRunCallBack}
callBackPath=${callBackPath}
"

runR="Rscript --no-save "

###############################################################################
#'                              Run Find Cluster                             '#
###############################################################################
if [ $toRunClustering = "YES" ]; then
    ${runR} ${HOME}/configs/public/pipeline/UMAP_CLUSTER_JOBS_EMBEDED/FindCluster.R -d ${inData} -o ${currentFolder}/cluster.rds -r ${reduction} -n ${npc} -e ${res}
fi



###############################################################################
#'                             Run Common Analysis                           '#
###############################################################################
srcD=${HOME}/configs/public/src
ResD=${currentFolder}
paramD=${HOME}/configs/public/params
databaseD=${HOME}/configs/public/knowledge/database

if [ $toRunCommonAnalysis = "YES" ]; then

    if [ ! -d ${ResD}/bubbleplot ]; then
        mkdir -p ${ResD}/bubbleplot
    fi

    ${runR} ${srcD}/bubble-plot.R -d ${ResD}/cluster.rds -o ${ResD}/bubbleplot -m ${databaseD}/TMarkers.txt
    ${runR} ${srcD}/bubble-plot.R -d ${ResD}/cluster.rds -o ${ResD}/bubbleplot -m ${databaseD}/topLevel.txt
    ${runR} ${srcD}/bubble-plot.R -d ${ResD}/cluster.rds -o ${ResD}/bubbleplot -m ${databaseD}/general/generalAll.txt

    if [ ! -d ${ResD}/featureplot ]; then
        mkdir -p ${ResD}/featureplot
    fi

    # ${runR} ${srcD}/feature-plot.R -d ${ResD}/cluster.rds -o ${ResD}/featureplot/topLevel -c ${paramD}/feature-plot-origin.json -m ${databaseD}/topLevel.txt
    # ${runR} ${srcD}/feature-plot.R -d ${ResD}/cluster.rds -o ${ResD}/featureplot/tmarkers -c ${paramD}/feature-plot-origin.json -m ${databaseD}/TMarkers.txt

    ${runR} ${srcD}/qc-by-cluster.R -d ${ResD}/cluster.rds -o ${ResD}/qc-by-cluster.pdf
    ${runR} ${srcD}/visualize.R -d ${ResD}/cluster.rds
    # ${runR} ${srcD}/visualize_batch.R -d ${ResD}/cluster.rds
    # ${runR} ${srcD}/findmarker.r -d ${ResD}/cluster.rds -o ${ResD}/snn-single-markers.tsv

    # ${runR} ${srcD}/snn-marker.R -d ${ResD}/cluster.rds -o ${ResD}/snn-markers.tsv -c ${paramD}/snn-marker.json
    # ${runR} ${srcD}/snn-heatmap.R -d ${ResD}/cluster.rds -o ${ResD}/markersHeatmap.pdf -c ${paramD}/snn-heatmap.json -m ${ResD}/snn-markers.tsv -p heatmap
    # python ${srcD}/statMarkers.py --markersTop ${ResD}/snn-markers.tsv --markersDatabase ${databaseD}/immuneCellMarkerAllinBox_Yanshuo.txt --out ${ResD}/markers.ys_celltype.tsv
    # python ${srcD}/statMarkers.py --markersTop ${ResD}/markers.top.tsv --markersDatabase ${databaseD}/immuneCellMarkerAllinBox_Yanshuo.txt --out ${ResD}/markers.top.ys_celltype.tsv
fi

###############################################################################
#'                             Run Callback Script                           '#
###############################################################################
if [ $toRunCallBack = "YES" ]; then
    if [ -f $callBackPath ]; then
        $callBackPath ${ResD}/cluster.rds
    fi
fi
