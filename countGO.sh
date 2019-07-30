#!/bin/bash

download_seqs()
{
    #Dataset name from Ensembl
    dataset=$1
    dataset=\"$dataset\"

    #GO terms must be separated by comma (no space)
    goterms=$2
    goterms=\"$goterms\"

    wget -o /dev/null -O $species.fa 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?> <!DOCTYPE Query> <Query virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" completionStamp = "1"> <Dataset name = '$dataset' interface = "default" > <Filter name = "go_parent_term" value = '$goterms'/> <Attribute name = "peptide" /> <Attribute name = "ensembl_peptide_id" /> <Attribute name = "ensembl_gene_id" /> <Attribute name = "ensembl_transcript_id" /> <Attribute name = "gene_biotype" /> <Attribute name = "external_gene_name" /> <Attribute name = "description" /> </Dataset> </Query>'
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

queryid=$2

input="data/dataSpecies.txt"
while read -r species; do
    echo "Reading data for "$species"."

    goterms=$1
    dataset=$(python get_dataset.py $species)
    download_seqs $dataset $goterms $species
    numTerms=$(grep -o ">" count/$species.fa | wc -l)
    echo "${numTerms}" >> count\_$queryid.txt
done < "$input"

rm -rf count

end=$SECONDS
elapsed=$(( end - start ))
convertsecs $elapsed