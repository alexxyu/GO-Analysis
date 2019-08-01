#!/bin/bash

download_genes()
{
    #Dataset name from Ensembl
    dataset=$1
    dataset=\"$dataset\"

    #GO terms must be separated by comma (no space)
    goterms=$2
    goterms=\"$goterms\"


    wget -o /dev/null -O $species.fa 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?> <!DOCTYPE Query> <Query virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6"> <Dataset name = '$dataset' interface = "default" > <Filter name = "go_parent_term" value = '$goterms'/> <Attribute name = "ensembl_gene_id" /> </Dataset> </Query>'
    mv $species.fa count/$species.fa
}

convertsecs()
{
    ((h=${1}/3600))
    ((m=(${1}%3600)/60))
    ((s=${1}%60))
    printf "Elapsed time for total process: %02d:%02d:%02d\n" $h $m $s
}

start=$SECONDS
mkdir -p count
mkdir -p countData

input="data/go_terms.txt"
species="Rat"

while read -r term; do

    dataset=$1
    if grep -q $term banList.txt; then
        continue
    else
    
        echo "Working on GO term" $term
        download_genes $dataset $term $species
        numTerms=$(wc -l count/$species.fa | awk '{print $1}')
        
        if [ "${numTerms% *}" -gt 120 ]; then
            #Remove ancestor terms since they will have more terms
            Rscript filter.R $term 1
        elif [ "${numTerms% *}" -lt 10 ]; then
            #Remove descendant terms since they will have fewer terms
            Rscript filter.R $term 0
        else
            echo "${term}" >> out.txt
        fi

        echo "${term}" >> banList.txt

    fi
    
done < "$input"

rm -rf count

end=$SECONDS
elapsed=$(( end - start ))
convertsecs $elapsed