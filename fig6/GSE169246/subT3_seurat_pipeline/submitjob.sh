#BSUB -J seurat_h1
#BSUB -q e80medium
#BSUB -w 'done(19238686)'
#BSUB -W 23:00
#BSUB -n 1
#BSUB -M 550
#BSUB -R rusage[mem=550]
#BSUB -B
#BSUB -N
#BSUB -u ychu2@mdanderson.org
#BSUB -o /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT3_seurat_pipeline/submitjob.o.txt
#BSUB -e /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT3_seurat_pipeline/submitjob.e.txt
#BSUB -cwd /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT3_seurat_pipeline/
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT3_seurat_pipeline/submitjob.o.txt
rm -rf /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT3_seurat_pipeline/submitjob.e.txt
module load python/3.7.3-anaconda
module load R/4.0.3
#____----____----____

/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/public/public/UMAP_CLUSTER_JOBS_EMBEDED/run.sh \
    --inData  /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/result/GSE169246/subT2_harmony/1_run_harmony/harmony_output.rds\
    --reduction harmony \
    --mainscriptsFolder  /rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT3_seurat_pipeline/ \
    --parentJobName "GSE169246_subT_h1" \
    --npcArray "15;20;30;50" \
    --UMAPDistArray "0.1" \
    --ClusterResArray "0.3;0.6;0.5" \
    --NneighborsArray "50" \
    --toRunUMAP "YES" \
    --toRunClustering "YES" \
    --toRunCommonAnalysis "YES" \
    --toRunCallBack "YES" \
    --callBackPath "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/R3Q7/code/pipeline/private/GSE169246/subT3_seurat_pipeline/callBack.sh"
