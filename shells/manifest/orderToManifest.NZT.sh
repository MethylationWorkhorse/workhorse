
PRG=orderToManifest.NZT

TOP=/Users/bbarnes/Documents/Projects/workhorse
EXE=${TOP}/scripts/manifest/${PRG}.pl

# AQP1:: NZT/orders/Design/selected.order.csv
#
NAM=NZT
ORD_1=${TOP}/dat/NZT/orders/Design/selected.order.csv
MAT_1=${TOP}/dat/NZT/AQP1/08-29-2018_11_32_24/20297484_probes.match.gz
AQP_1=${TOP}/dat/NZT/AQP1/BS0031918-AQP.txt.gz

# AQP2::
#
NAM=NZT
ORD_2=${TOP}/dat/NZT/orders/NZT-round2.orders.FINAL.Jan2019.csv.gz
MAT_2=${TOP}/dat/NZT/AQP2/20297484_probes.match.gz
AQP_2=${TOP}/dat/NZT/AQP2/BS0032272-AQP.txt.gz

BAM=${TOP}/dat/NZT/alignment/nzt.order.bsmap.hg37.sorted.bam
OUT=${TOP}/workspace/${PRG}/org/${NAM}

BAM=${TOP}/aln/NZT-BP2.manifest.seq.genomic.sorted.bam
OUT=${TOP}/workspace/${PRG}/bsp/${NAM}

mkdir -p ${OUT}

echo "EXE="${EXE}
echo "OUT="${OUT}

# ${EXE} -v 5 -srd \
${EXE} -v 5 \
    -out ${OUT} \
    -nam ${NAM} \
    -ord ${ORD_1} -mat ${MAT_1} -aqp ${AQP_1} \
    -ord ${ORD_2} -mat ${MAT_2} -aqp ${AQP_2} \
    -bam ${BAM}


# End of file
