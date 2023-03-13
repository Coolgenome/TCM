#BSUB -J p1CD4_V7
#BSUB -q short
#BSUB -W 1:00
#BSUB -n 1
#BSUB -M 10
#BSUB -R rusage[mem=10]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/CD4/p1CD4_V7.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/CD4/p1CD4_V7.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/CD4/
rm -rf /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/CD4/p1CD4_V7.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/CD4/p1CD4_V7.e.txt
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


ResD=${analysisPath}/validate/CD4_V7
if [ ! -d $ResD ]; then
    mkdir -p $ResD
fi

runR="Rscript --no-save "


##do the job
echo "load data"

#' submit multile jobs to normalize and scale separately #############################################
# ${runR} ${srcD}/load-scaleTObjectList.R -d ${DataD} -o ${ResD}/PCA.rds

# ${runR} ${srcD}/RunPCA_RPCA.R -d ${ResD}/data.rds -o ${ResD}/pca.rds
# ${runR} ${srcD}/visualize_PCAgenes.R -d ${ResD}/pca.rds

${HOME}/share/UMAP_CLUSTER_JOBS_EMBEDED/run.sh \
       --inData ${ResD}/pca.rds \
       --reduction pca \
       --mainscriptsFolder ${mainscriptsFolder} \
       --parentJobName "p1CD4_V7" \
       --npcArray "50" \
       --UMAPDistArray "0.1" \
       --ClusterResArray "0.3" \
       --NneighborsArray "50" \
       --toRunUMAP "NO" \
       --toRunClustering "NO" \
       --toRunCommonAnalysis "NO" \
       --toRunCallBack "YES" \
       --callBackPath "/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/CD4/callBack_CD4.sh"

