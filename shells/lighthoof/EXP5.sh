
verbosity=4

EXP_NUM=EXP5

fresh=false
cluster=true
parallel=true
autoDetect=true
single=false

writeSigs=false
writeCall=false
writeAuto=true
plotAuto=false
minPval=0.02

platform=EPIC
manifest=B4
prgmTag=workhorse
EXE_NAME=lighthoof_main

# TOP_MAC=/Users/bbarnes/Documents/Projects/${prgmTag}
TOP_MAC=/Users/bbarnes/Documents/CustomerFacing
TOP_LIX=/illumina/scratch/darkmatter/data

if [ -e ${TOP_MAC} ]; then
    TOP=${TOP_MAC}
    DAT=${TOP_MAC}/dat
    SRC=${TOP_MAC}/git/${prgmTag}
    CONDA=mac
    RSCRIPT=Rscript
elif [ -e ${TOP_LIX} ]; then
    TOP=${TOP_LIX}
    DAT=${TOP_LIX}/dat
    SRC=${TOP_LIX}/git/${prgmTag}
    CONDA=conda_4.6.8
    # CONDA=Anaconda2-2019.10-Linux-x86_64
    # CONDA=Anaconda3-2019.10-Linux-x86_64
    RSCRIPT=/illumina/scratch/darkmatter/thirdparty/${CONDA}/bin/Rscript
else
    echo "Unrecognized top directory!"
    exit
fi

EXE=${SRC}/scripts/R/${EXE_NAME}.R
outDir=${TOP}/workspace/${EXE_NAME}/builds/${EXP_NUM}
idatsDir=${TOP}/idats_${EXP_NUM}
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
