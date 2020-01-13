
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 idat_name"
    exit 1
fi

# idatName=idats_ref
# idatName=DeltaBetaCore
idatName=$1

prgmTag=workhorse

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

# Commonly Modified Parameters::
single=false
parallel=true
cluster=true
verbosity=3

R_DIR=${SRC}/scripts/R
EXE=${R_DIR}/${prgmTag}.R

# Directory Parameters::
outDir=${TOP}/builds/${prgmTag}/${idatName}
idatsDir=${TOP}/idats_${idatName}
srcDir=${SRC}/scripts
funcScript=${srcDir}/R/workhorse_functions.R

# Generally Default Parameter
autoDetect=true
writeSSheet=true
writeMIDMAN=true
writeMAN=true
writeRDS=true
writeCSV=false
loadMAN=true
loadRDS=true
saveRDS=true
overRDS=true

platform=EPIC-B4
build=hg19
method=both
poobMinPval=0.2
negsMinPval=0.02

# Calculated internally by script::
#    --prgmTag=${prgmTag} \
#    --datDir=${datDir} \
#    --topDir=${topDir} \
#    --runMode=${runMode} \

# Debug Defaults (only use in interactive mode)::
#    --retSSET=false \
#    --retPRBS=false \

# echo "Rscript="${RSCRIPT}
# echo "exe="${EXE}

CMD=${RSCRIPT}" "${EXE}
CMD+=" --"verbosity=${verbosity}
CMD+=" --"Rscript=${RSCRIPT}
CMD+=" --"topDir=${TOP}
CMD+=" --"outDir=${outDir}
CMD+=" --"idatsDir=${idatsDir}
CMD+=" --"srcDir=${srcDir}
CMD+=" --"funcScript=${funcScript}
CMD+=" --"platform=${platform}
CMD+=" --"build=${build}
CMD+=" --"method=${method}
CMD+=" --"poobMinPval=${poobMinPval}
CMD+=" --"negsMinPval=${negsMinPval}
if [ "${autoDetect}" = true ]; then
    CMD+=" --autoDetect"
fi
if [ "${writeSSheet}" = true ]; then
    CMD+=" --writeSSheet"
fi
if [ "${writeMIDMAN}" = true ]; then
    CMD+=" --writeMIDMAN"
fi
if [ "${writeMAN}" = true ]; then
    CMD+=" --writeMAN"
fi
if [ "${writeRDS}" = true ]; then
    CMD+=" --writeRDS"
fi
if [ "${writeCSV}" = true ]; then
    CMD+=" --writeCSV"
fi
if [ "${loadMAN}" = true ]; then
    CMD+=" --loadMAN"
fi
if [ "${loadRDS}" = true ]; then
    CMD+=" --loadRDS"
fi
if [ "${saveRDS}" = true ]; then
    CMD+=" --saveRDS"
fi
if [ "${overRDS}" = true ]; then
    CMD+=" --overRDS"
fi
if [ "${single}" = true ]; then
    CMD+=" --single"
fi
if [ "${parallel}" = true ]; then
    CMD+=" --parallel"
fi
if [ "${cluster}" = true ]; then
    CMD+=" --cluster"
fi

mkdir -p ${outDir}

echo ${CMD}

${CMD}

## End of file
