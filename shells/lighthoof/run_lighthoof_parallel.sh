#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 idat_name"
    exit 1
fi
idatsDir=$1
outSuffix=$2
EXP_NAME=`basename $idatsDir`
EXP_NAME=$(sed 's/^idats_//' <<< "$EXP_NAME")

verbosity=4

fresh=false
cluster=true
parallel=true
autoDetect=true
single=false

buildSubDir=true
addBeadCounts=true

DyeSwapNoob=true
SwapOpen=true
RawOpen=false

writeSigs=true
writeCall=true
writeAuto=false
plotAuto=false
minPval=0.02

platform=EPIC
manifest=B4
prgmTag=workhorse
EXE_NAME=lighthoof_main

TOP_MAC=/Users/bbarnes/Documents/CustomerFacing
TOP_LIX=/illumina/scratch/darkmatter/data

if [ -e ${TOP_MAC} ]; then
    TOP=${TOP_MAC}
    SRC=${TOP_MAC}/${prgmTag}
    DAT=${SRC}/dat
    CONDA=mac
    RSCRIPT=/usr/local/bin/Rscript
elif [ -e ${TOP_LIX} ]; then
    TOP=${TOP_LIX}
    SRC=${TOP_LIX}/${prgmTag}
    DAT=${SRC}/dat
    CONDA=conda_4.6.8
    # CONDA=Anaconda2-2019.10-Linux-x86_64
    # CONDA=Anaconda3-2019.10-Linux-x86_64
    RSCRIPT=/illumina/scratch/darkmatter/thirdparty/${CONDA}/bin/Rscript
else
    echo "Unrecognized top directory!"
    exit
fi

EXE=${SRC}/scripts/R/${EXE_NAME}.R
outDir=${TOP}/workspace/${EXE_NAME}${outSuffix}_par/builds/${EXP_NAME}
datDir=${TOP}/dat

CMD=${RSCRIPT}" "${EXE}
CMD+=" --"Rscript=${RSCRIPT}
CMD+=" --"outDir=${outDir}
CMD+=" --"datDir=${datDir}
CMD+=" --"idatsDir=${idatsDir}
CMD+=" --"platform=${platform}
CMD+=" --"manifest=${manifest}
CMD+=" --"minPval=${minPval}
if [ "${fresh}" = true ]; then
    CMD+=" --fresh"
fi
if [ "${single}" = true ]; then
    CMD+=" --single"
fi
if [ "${cluster}" = true ]; then
    CMD+=" --cluster"
fi
if [ "${parallel}" = true ]; then
    CMD+=" --parallel"
fi
if [ "${autoDetect}" = true ]; then
    CMD+=" --autoDetect"
fi

if [ "${buildSubDir}" = true ]; then
    CMD+=" --buildSubDir"
fi
if [ "${addBeadCounts}" = true ]; then
    CMD+=" --addBeadCounts"
fi

if [ "${DyeSwapNoob}" = true ]; then
    CMD+=" --DyeSwapNoob"
fi
if [ "${SwapOpen}" = true ]; then
    CMD+=" --SwapOpen"
fi
if [ "${RawOpen}" = true ]; then
    CMD+=" --RawOpen"
fi


if [ "${writeSigs}" = true ]; then
    CMD+=" --writeSigs"
fi
if [ "${writeCall}" = true ]; then
    CMD+=" --writeCall"
fi
if [ "${writeAuto}" = true ]; then
    CMD+=" --writeAuto"
fi
if [ "${plotAuto}" = true ]; then
    CMD+=" --plotAuto"
fi
CMD+=" --"verbosity=${verbosity}

# mkdir -p ${outDir}
echo ${CMD}

${CMD}

## End of file
