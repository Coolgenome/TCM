##!/usr/bin/env bash
#' filename : zoomIn.sh
#' Date : 2022-08-15
#' contributor : Yanshuo Chu
#' function: convert cut images and add a map scale

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        --inDataPath)
            inDataPath="$2"
            shift # past argument
            shift # past value
            ;;
        --outDataPath)
            outDataPath="$2"
            shift # past argument
            shift # past value
            ;;
        --xmin)
            xmin="$2"
            shift # past argument
            shift # past value
            ;;
        --xmax)
            xmax="$2"
            shift # past argument
            shift # past value
            ;;
        --ymin)
            ymin="$2"
            shift # past argument
            shift # past value
            ;;
        --ymax)
            ymax="$2"
            IFS=';' read -a npcArray <<< "$npcArray"
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
inDataPath=${inDataPath}
outDataPath=${outDataPath}
xmin=${xmin}
xmax=${xmax}
ymin=${ymin}
ymax=${ymax}
"

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 "$1"
fi

convert $inDataPath -crop $((${xmax}-${xmin}))x$((${ymax}-${ymin}))+${xmin}+${ymin} ${outDataPath}.bak

module load R/3.6.0
echo " Rscript --no-save /rsrch3/scratch/genomic_med/ychu2/projects/p1review/CosMx/code/pipeline/private/1_MapScale/plotMapScale.R -d ${outDataPath}.bak -o ${outDataPath} -w $((${xmax}-${xmin})) -h $((${ymax}-${ymin})) "
Rscript --no-save /rsrch3/scratch/genomic_med/ychu2/projects/p1review/CosMx/code/pipeline/private/1_MapScale/plotMapScale.R -d ${outDataPath}.bak -o ${outDataPath} -w $((${xmax}-${xmin})) -h $((${ymax}-${ymin}))

module load R/4.1.0
echo " Rscript --no-save /rsrch3/scratch/genomic_med/ychu2/projects/p1review/CosMx/code/pipeline/private/1_MapScale/plotMapScale.R -d ${outDataPath}.bak -o ${outDataPath} -w $((${xmax}-${xmin})) -h $((${ymax}-${ymin})) "
Rscript --no-save /rsrch3/scratch/genomic_med/ychu2/projects/p1review/CosMx/code/pipeline/private/1_MapScale/plotGiotto.R \
        -f $(echo $(basename $inDataPath) | egrep -o "[0-9]+") --xmin ${xmin} --xmax ${xmax} --ymin ${ymin} --ymax ${ymax} \
        -o ${outDataPath%.jpg}_loc.png


