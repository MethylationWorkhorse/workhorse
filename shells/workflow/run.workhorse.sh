
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 idat_name"
    exit 1
fi

# idatName=idats_ref
# idatName=DeltaBetaCore
idatName=$1

prgmTag=workhorse

# TOP_MAC=/Users/bbarnes/Documents/Projects/${prgmTag}
TOP_MAC=/Users/bbarnes/Documents/Projects/${prgmTag}/git/workhorse-master
TOP_LIX=/illumina/scratch/darkmatter/Projects/${prgmTag}

if [ -e ${TOP_MAC} ]; then
    TOP_SRC=${TOP_MAC}
    TOP_DIR=/Users/bbarnes/Documents/Projects/${prgmTag}
    CONDA=mac
    Rscript=/usr/local/bin/Rscript
elif [ -e ${TOP_LIX} ]; then
    TOP_SRC=${TOP_LIX}
    TOP_DIR=${TOP_LIX}
    CONDA=conda_4.6.8
    # CONDA=Anaconda2-2019.10-Linux-x86_64
    # CONDA=Anaconda3-2019.10-Linux-x86_64
    Rscript=/illumina/scratch/darkmatter/thirdparty/${CONDA}/bin/Rscript
else
    echo "Unrecognized top directory!"
    exit
fi

# Commonly Modified Parameters::
single=false
parallel=true
cluster=true
verbosity=3

# Source Directory Parameters::
srcDir=${TOP_SRC}/scripts
funcScript=${srcDir}/R/workhorse_functions.R
EXE=${srcDir}/R/${prgmTag}.R

# Data Directory Parameters::
topDir=${TOP_DIR}
outDir=${TOP_DIR}/builds/${prgmTag}/${idatName}
idatsDir=${TOP_DIR}/idats/${idatName}

# Generally Default Parameter
autoDetect=true
writeSSheet=true
writeMIDMAN=true
writeMAN=true
writeGRP=true
writeRDS=true
writeCSV=false
loadMAN=false
loadGRP=false
loadRDS=true
saveRDS=true
overRDS=true

if [ 1 ]; then
    autoDetect=true
    writeSSheet=true
    writeMIDMAN=false
    writeMAN=false
    writeGRP=false
    writeRDS=true
    writeCSV=false
    loadMAN=true
    loadGRP=true
    loadRDS=true
    saveRDS=true
    overRDS=true

    echo "Low Ouput"
fi

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

CMD=${Rscript}" "${EXE}
CMD+=" --"verbosity=${verbosity}
CMD+=" --"Rscript=${Rscript}
CMD+=" --"outDir=${outDir}
CMD+=" --"idatsDir=${idatsDir}
CMD+=" --"topDir=${topDir}
# CMD+=" --"srcDir=${srcDir}
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
if [ "${writeGRP}" = true ]; then
    CMD+=" --writeGRP"
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
if [ "${loadGRP}" = true ]; then
    CMD+=" --loadGRP"
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
