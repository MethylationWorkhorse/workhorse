

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 name shell"
    exit 1
fi
qsub -cwd -pe threaded 16 -l excl=true -N $1 $2
