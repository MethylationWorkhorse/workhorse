

TOP=/Users/bbarnes/Documents/Projects/workhorse
EXE=${TOP}/scripts/copy/copy-check.pl

IDATS=/Users/bbarnes/Documents/Projects/darkmatter/Projects/sesamize/idats_DeltaBetaCore

NAM=DeltaBetaCore
OUT=${TOP}/idats/${NAM}

ls ${IDATS}/*/*.idat.gz | ${EXE} -o ${OUT}

# End of file
