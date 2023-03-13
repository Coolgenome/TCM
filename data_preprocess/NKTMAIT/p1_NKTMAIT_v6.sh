#BSUB -J p1_NKTMAIT_v6
#BSUB -q highmem
#BSUB -W 24:00
#BSUB -n 1
#BSUB -M 500
#BSUB -R rusage[mem=500]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/NKTMAIT/p1_NKTMAIT_v6.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/NKTMAIT/p1_NKTMAIT_v6.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/NKTMAIT/
rm -rf /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/NKTMAIT/p1_NKTMAIT_v6.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/NKTMAIT/p1_NKTMAIT_v6.e.txt
module load python/3.7.3-anaconda
module load R/3.6.0
#____----____----____----____----____----____----____----____----____----____----____----____----____----____----____----____----



projectPath=/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject
DataD=${projectPath}/data/T/filterd

analysisPath=${projectPath}/analysis
mainscriptsFolder=${analysisPath}/scripts
pipelinesFolder=${mainscriptsFolder}/pipelines
srcD=${analysisPath}/scripts/src
paramD=${analysisPath}/scripts/params
databaseD='/rsrch3/home/genomic_med/ychu2/share/database'

ResD=${analysisPath}/validate/NKTMAIT_V6
if [ ! -d $ResD ]; then
    mkdir -p $ResD
fi

runR="Rscript --no-save "


##do the job
echo "load data"

#' submit multile jobs to normalize and scale separately #############################################
#${runR} ${srcD}/load-scaleTObjectList.R -d ${DataD} -o ${ResD}/PCA.rds

###############################################################################
#'                              split subcluster                             '#
# ${runR} ${srcD}/split-cluster.R -i /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/validate/CD8_V6/nPC_50/UMAP_dist_0.1_nneighbor_50/p1CD8_V6_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds -o ${ResD}/data.rds -c "7;1"
###############################################################################

# ${runR} ${srcD}/RunPCA_RPCA.R -d ${ResD}/data.rds -o ${ResD}/pca.rds
# ${runR} ${srcD}/visualize_PCAgenes.R -d ${ResD}/pca.rds

${HOME}/share/UMAP_CLUSTER_JOBS_EMBEDED/run.sh \
       --inData ${ResD}/pca.rds \
       --reduction pca \
       --mainscriptsFolder ${mainscriptsFolder} \
       --parentJobName "p1_NKTMAIT_v6" \
       --npcArray "5" \
       --UMAPDistArray "0.1" \
       --ClusterResArray "0.3" \
       --NneighborsArray "35" \
       --toRunUMAP "NO" \
       --toRunClustering "NO" \
       --toRunCommonAnalysis "NO" \
       --toRunCallBack "YES" \
       --callBackPath "${pipelinesFolder}/split/NKTMAIT/callBack_NKTMAIT.sh"

