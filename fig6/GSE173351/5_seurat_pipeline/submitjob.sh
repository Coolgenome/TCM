#BSUB -J submitjob
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 1
#BSUB -M 10
#BSUB -R rusage[mem=10]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/5_seurat_pipeline/submitjob.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/5_seurat_pipeline/submitjob.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/5_seurat_pipeline/
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/5_seurat_pipeline/submitjob.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/5_seurat_pipeline/submitjob.e.txt
module load python/3.7.3-anaconda
module load R/3.6.0
#____----____----____


/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/public/public/UMAP_CLUSTER_JOBS_EMBEDED/run.sh \
    --inData /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE173351/4_harmony/harmony_output.rds \
    --reduction harmony \
    --mainscriptsFolder  /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/5_seurat_pipeline \
    --parentJobName "GSE173351" \
    --npcArray "30" \
    --UMAPDistArray "0.1" \
    --ClusterResArray "0.3" \
    --NneighborsArray "50" \
    --toRunUMAP "NO" \
    --toRunClustering "NO" \
    --toRunCommonAnalysis "YES" \
    --toRunCallBack "NO" \
    --callBackPath "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE173351/5_seurat_pipeline/callBack.sh"
