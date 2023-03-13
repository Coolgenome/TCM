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

# ${runR} ${srcD}/iml_check_bubble-plot.R -d $dataPath -c "NKT_c0_CD8 MAIT-like_c1 MAIT_c2 MAIT-like_c3_SELL NKT_c4_KIR2D"

${runR} ${srcD}/bubble-plot.R -d $dataPath -o $(dirname ${dataPath})/bubbleplot -m  /rsrch3/home/genomic_med/ychu2/share/database/Pan-T/Reference/PMID34290406.txt

# ${runR} ${srcD}/AlluvialPlot.R -l /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/NKGDT_V5/nPC_10/UMAP_dist_0.1_nneighbor_35/p1NKGDT_V4_10_UMAP_dist_0.1_nneighbor_35_CLUSTER_res_0.3/cluster.rds -r $dataPath
# ${runR} ${srcD}/DEG-bubble-plot.R -d $dataPath --DEGs $(dirname ${dataPath})/snn-single-markers.tsv --number 200
# ${runR} ${srcD}/tissue-composition-plot.R -d $dataPath


# FracArray=($(seq 0.1 0.1 0.8))
# for tempFrac in "${FracArray[@]}"; do
#     JOBNAME=cd8_c1c7_monocle2_${tempFrac}
#     tempFolder=$(dirname $dataPath)
#     JOBFOLDER=${tempFolder}
#     if [ ! -d $tempFolder ]; then
#         mkdir -p $tempFolder
#     fi
#     if [ -f ${JOBFOLDER}/${JOBNAME}.o.txt ] || [ -f ${JOBFOLDER}/${JOBNAME}.e.txt ]; then
#         rm ${JOBFOLDER}/${JOBNAME}.*.txt -f
#     fi
#     bsub \
#         -J ${JOBNAME} \
#         -o ${JOBFOLDER}/${JOBNAME}.o.txt \
#         -e ${JOBFOLDER}/${JOBNAME}.e.txt \
#         -cwd ${JOBFOLDER} \
#         -q long \
#         -W 120:00 \
#         -n 1 \
#         -M 800\
#         -R rusage[mem=800] \
#         -B \
#         -N \
#         -u ychu2@mdanderson.org \
#         /bin/bash -c "module load python/3.7.3-anaconda; module load R/3.6.0; ${runR} ${srcD}/monocle2.R -d ${dataPath} -f ${tempFrac}"
# done

# JOBNAME=cd8_c1c7_monocle3
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
#     /bin/bash -c "module load python/3.7.3-anaconda; module load R/3.6.0; ${runR} ${srcD}/trajectory_monocle3_from_seurat.r -d ${dataPath} -s $(dirname ${dataPath})/startCells.txt"
