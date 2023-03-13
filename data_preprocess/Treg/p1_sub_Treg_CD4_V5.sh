#BSUB -J p1_sub_Treg_CD4_V5
#BSUB -q short
#BSUB -W 1:00
#BSUB -n 1
#BSUB -M 1
#BSUB -R rusage[mem=1]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/p1_sub_Treg_CD4_V5.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/p1_sub_Treg_CD4_V5.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/
rm -rf /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/p1_sub_Treg_CD4_V5.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/scripts/pipelines/split/p1_sub_Treg_CD4_V5.e.txt
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

ResD=${analysisPath}/validate/Sub_Treg_CD4_V5
if [ ! -d $ResD ]; then
    mkdir -p $ResD
fi

runR="Rscript --no-save "


##do the job
echo "load data"

#' submit multile jobs to normalize and scale separately #############################################
#${runR} ${srcD}/load-scaleTObjectList.R -d ${DataD} -o ${ResD}/PCA.rds

# ${runR} ${srcD}/RunPCA_RPCA.R -d ${ResD}/data.rds -o ${ResD}/pca.rds
# ${runR} ${srcD}/visualize_PCAgenes.R -d ${ResD}/pca.rds

${HOME}/share/UMAP_CLUSTER_JOBS_EMBEDED/run.sh \
       --inData ${ResD}/pca.rds \
       --reduction pca \
       --mainscriptsFolder ${mainscriptsFolder} \
       --parentJobName "p1_sub_Treg_CD4_V5" \
       --npcArray "15" \
       --UMAPDistArray "0.1" \
       --ClusterResArray "0.3" \
       --NneighborsArray "20" \
       --toRunUMAP "NO" \
       --toRunClustering "NO" \
       --toRunCommonAnalysis "YES" \
       --toRunCallBack "YES" \
       --callBackPath "${pipelinesFolder}/split/callBack_Treg.sh"

       # --callBackPath "${pipelinesFolder}/split/callBack.sh"

