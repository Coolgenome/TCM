##!/usr/bin/env bash

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        --inData)
            inData="$2"
            shift # past argument
            shift # past value
            ;;
        --reduction)
            reduction="$2"
            shift # past argument
            shift # past value
            ;;
        --mainscriptsFolder)
            mainscriptsFolder="$2"
            shift # past argument
            shift # past value
            ;;
        --parentJobName)
            parentJobName="$2"
            shift # past argument
            shift # past value
            ;;
        --npcArray)
            npcArray="$2"
            IFS=';' read -a npcArray <<< "$npcArray"
            shift # past argument
            shift # past value
            ;;
        --UMAPDistArray)
            UMAPDistArray="$2"
            IFS=';' read -a UMAPDistArray <<< "$UMAPDistArray"
            shift # past argument
            shift # past value
            ;;
        --ClusterResArray)
            ClusterResArray="$2"
            shift # past argument
            shift # past value
            ;;
        --NneighborsArray)
            NneighborsArray="$2"
            IFS=';' read -a NneighborsArray <<< "$NneighborsArray"
            shift # past argument
            shift # past value
            ;;
        --toRunUMAP)
            toRunUMAP="$2"
            shift # past argument
            shift # past value
            ;;
        --toRunClustering)
            toRunClustering="$2"
            shift # past argument
            shift # past value
            ;;
        --toRunCommonAnalysis)
            toRunCommonAnalysis="$2"
            shift # past argument
            shift # past value
            ;;
        --toRunCallBack)
            toRunCallBack="$2"
            shift # past argument
            shift # past value
            ;;
        --callBackPath)
            callBackPath="$2"
            shift # past argument
            shift # past value
            ;;
        *)    # unknown option
            POSITIONAL+=("$1") # save it in an array for later
            shift # past argument
            ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo "run parameters:
inData=${inData}
reduction=${reduction}
mainscriptsFolder=${mainscriptsFolder}
parentJobName=${parentJobName}
npcArray[@]=${npcArray[@]}
UMAPDistArray[@]=${UMAPDistArray[@]}
ClusterResArray=${ClusterResArray}
NneighborsArray[@]=${NneighborsArray[@]}
toRunUMAP=${toRunUMAP}
toRunClustering=${toRunClustering}
toRunCommonAnalysis=${toRunCommonAnalysis}
toRunCallBack=${toRunCallBack}
callBackPath=${callBackPath}
"

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 "$1"
fi

rootFolder=$(dirname ${inData})

mmForInData=300
if [ "$mmForInData" -gt "1200" ]; then
    mmForInData=350
fi

qName="highmem"
wTime="2:00"
if [ "$mmForInData" -lt "80" ]; then
    qName="e80short"
    wTime="2:00"
else if [ "$mmForInData" -lt "100" ]; then
         qName="e80medium"
         wTime="4:00"
     else if [ "$mmForInData" -lt "200" ]; then
              qName="highmem"
              wTime="6:00"
          else if [ "$mmForInData" -lt "500" ]; then
                   qName="highmem"
                   wTime="12:00"
               else if [ "$mmForInData" -lt "1000" ]; then
                        qName="highmem"
                        wTime="48:00"
                    else if [ "$mmForInData" -lt "2000" ]; then
                             qName="highmem"
                             wTime="120:00"
                         fi
                    fi
               fi
          fi
     fi
fi


if [ "$toRunUMAP" == "NO" ]; then
    mmForInData=50
    qName="e80short"
    wTime="2:50"
fi

for npc in "${npcArray[@]}"; do
    npcFolder=${rootFolder}/nPC_${npc}
    for dist in "${UMAPDistArray[@]}"; do
        for nneighbor in "${NneighborsArray[@]}"; do
            tempFolder=${npcFolder}/UMAP_dist_${dist}_nneighbor_${nneighbor}
            if [ ! -d $tempFolder ]; then
                mkdir -p $tempFolder
            fi
            JOBFOLDER=${tempFolder}
            JOBNAME=${parentJobName}_UMAP_dist_${dist}_nneighbor_${nneighbor}
            if [ -f ${JOBFOLDER}/${JOBNAME}.o.txt ] || [ -f ${JOBFOLDER}/${JOBNAME}.e.txt ]; then
                rm ${JOBFOLDER}/${JOBNAME}.*.txt -f
            fi

            bsub \
                -J ${JOBNAME} \
                -o ${JOBFOLDER}/${JOBNAME}.o.txt \
                -e ${JOBFOLDER}/${JOBNAME}.e.txt \
                -cwd ${JOBFOLDER} \
                -q ${qName} \
                -W ${wTime} \
                -n 1 \
                -M ${mmForInData} \
                -R rusage[mem=${mmForInData}] \
                -B \
                -N \
                -u ychu2@mdanderson.org \
                /bin/bash -c "${HOME}/configs/public/pipeline/UMAP_CLUSTER_JOBS_EMBEDED/RunUMAPJobs.sh ${mainscriptsFolder} ${inData} ${tempFolder} ${dist} ${reduction} ${npc} ${JOBNAME} \"${ClusterResArray}\" ${nneighbor} ${toRunUMAP} ${toRunClustering} ${toRunCommonAnalysis} ${toRunCallBack} ${callBackPath}"
        done
    done
done
