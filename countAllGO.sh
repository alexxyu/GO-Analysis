#!/bin/bash

input="out.txt"
opSize=$(wc -l out.txt | awk '{print $1}')
counter=1
while read -r term; do
    echo "GO term number" $counter "out of" $opSize"."
    echo "Working on GO term" $term"."; echo;
    bash countGO.sh $term ${term:3}
done < "$input"