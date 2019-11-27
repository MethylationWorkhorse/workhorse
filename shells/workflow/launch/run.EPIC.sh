
idatName=EPIC
prgmTag=workhorse

TOP_MAC=/Users/bbarnes/Documents/Projects/${prgmTag}/git/${prgmTag}
TOP_LIX=/illumina/scratch/darkmatter/Projects/${prgmTag}/git/${prgmTag}

if [ -e ${TOP_MAC} ]; then
    TOP=${TOP_MAC}
    CONDA=mac
    RSCRIPT=Rscript
elif [ -e ${TOP_LIX} ]; then
    TOP=${TOP_LIX}
    CONDA=conda_4.6.8
    # CONDA=Anaconda2-2019.10-Linux-x86_64
    # CONDA=Anaconda3-2019.10-Linux-x86_64
    RSCRIPT=/illumina/scratch/darkmatter/thirdparty/${CONDA}/bin/Rscript
else
    echo "Unrecognized top directory!"
    exit
fi

EXE=${TOP}/shells/workflow/run.workhorse.sh

${EXE} ${idatName}

## End of file
