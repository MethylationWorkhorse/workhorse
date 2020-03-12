#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 idat_name"
    exit 1
fi
idatsDir=$1
outSuffix=$2
EXP_NAME=`basename $idatsDir`
EXP_NAME=$(sed 's/^idats_//' <<< "$EXP_NAME")

prgmTag=workhorse
EXE_NAME=swifthoof_main

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

CMD=${RSCRIPT}" "${EXE}
CMD+=" --"Rscript=${RSCRIPT}

# Directories::
outDir=${TOP}/workspace/${EXE_NAME}${outSuffix}_par/builds/${EXP_NAME}
datDir=${TOP}/${prgmTag}/dat
auto_sam_csv=${datDir}/ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_DELTA-Zymo.csv.gz
# idatsDir=NULL
CMD+=" --auto_sam_csv"=${auto_sam_csv}
CMD+=" --outDir"=${outDir}
CMD+=" --datDir"=${datDir}
CMD+=" --idatsDir"=${idatsDir}

# Optional Files::
subManifest=false
# manifestPath=NULL
# addressPath=NULL
if [ "${subManifest}" = true ]; then
    CMD+=" --subManifest"
fi
# CMD+=" --manifestPath"=${manifestPath}
# CMD+=" --addressPath"=${addressPath}

# Platform/Method Options::
platform='EPIC'
manifest='B4'
CMD+=" --platform"=${platform}
CMD+=" --manifest"=${manifest}

# Run Options::
fresh=false
buildSubDir=true
autoDetect=true
workflows='ind'
if [ "${fresh}" = true ]; then
    CMD+=" --fresh"
fi
if [ "${buildSubDir}" = true ]; then
    CMD+=" --buildSubDir"
fi
if [ "${autoDetect}" = true ]; then
    CMD+=" --autoDetect"
fi
CMD+=" --workflows"=${workflows}

# Output Options::
loadIdat=true
saveIdat=true
if [ "${loadIdat}" = true ]; then
    CMD+=" --loadIdat"
fi
if [ "${saveIdat}" = true ]; then
    CMD+=" --saveIdat"
fi

loadSsets=true
saveSsets=true
saveRawSset=false
if [ "${loadSsets}" = true ]; then
    CMD+=" --loadSsets"
fi
if [ "${saveSsets}" = true ]; then
    CMD+=" --saveSsets"
fi
if [ "${saveRawSset}" = true ]; then
    CMD+=" --saveRawSset"
fi

writeSset=false
writeCalls=true
writeSsheet=true
writeAuto=false
if [ "${writeSset}" = true ]; then
    CMD+=" --writeSset"
fi
if [ "${writeCalls}" = true ]; then
    CMD+=" --writeCalls"
fi
if [ "${writeSsheet}" = true ]; then
    CMD+=" --writeSsheet"
fi
if [ "${writeAuto}" = true ]; then
    CMD+=" --writeAuto"
fi

# Reporting Options::
sigs_sum_field='avg'
CMD+=" --sigs_sum_field"=${sigs_sum_field}

# Threshold Options::
minNegPval=0.02
minOobPval=0.1
minNegPerc=96
minOobPerc=85
minDeltaBeta=0.2
CMD+=" --minNegPval"=${minNegPval}
CMD+=" --minOobPval"=${minOobPval}
CMD+=" --minNegPerc"=${minNegPerc}
CMD+=" --minOobPerc"=${minOobPerc}
CMD+=" --minDeltaBeta"=${minDeltaBeta}

percisionBeta=4
percisionPval=6
CMD+=" --percisionBeta"=${percisionBeta}
CMD+=" --percisionPval"=${percisionPval}

# Parallel/Cluster Options::
single=false
parallel=true
cluster=true
if [ "${single}" = true ]; then
    CMD+=" --single"
fi
if [ "${parallel}" = true ]; then
    CMD+=" --parallel"
fi
if [ "${cluster}" = true ]; then
    CMD+=" --cluster"
fi

# Plotting Options::
plotSset=false
plotCalls=false
plotAuto=false
if [ "${plotSset}" = true ]; then
    CMD+=" --plotSset"
fi
if [ "${plotCalls}" = true ]; then
    CMD+=" --plotCalls"
fi
if [ "${plotAuto}" = true ]; then
    CMD+=" --plotAuto"
fi

plotFormat='pdf'
plotFormat='png'
CMD+=" --plotFormat"=${plotFormat}

dpi=72
dpi=120
CMD+=" --dpi"=${dpi}

plotMax=10000
plotSub=5000
CMD+=" --plotMax"=${plotMax}
CMD+=" --plotSub"=${plotSub}

# Verbosity Options::
verbosity=3
CMD+=" --"verbosity=${verbosity}

# mkdir -p ${outDir}
echo ${CMD}

${CMD}

## End of file
