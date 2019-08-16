#!/bin/bash

#Must have wget installed.

if [ $# -lt 2 ]; then
    echo "Error -- Program takes arguments as follows: "
    echo "bash countAllGO.sh [Path to file containing GO terms] [output filename]"
    exit 1
fi

input=$1
opSize=$(wc -l $input | awk '{print $1}')
counter=1

while read -r term; do
    echo "Working on GO term" $term "("$counter "of" $opSize")."
    bash countGO.sh $term ${term:3}
    let counter++
    echo
done < "$input"

mkdir -p output
python GOAnalysis.py $2
