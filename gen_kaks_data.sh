#!/bin/bash

download_seqs()
{
    #Dataset name from Ensembl
    dataset=$1
    dataset=\"$dataset\"

    #GO terms must be separated by comma (no space)
    goterms=$2
    goterms=\"$goterms\"

    filename=$3
    wget -O $3 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?> <!DOCTYPE Query> <Query virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" completionStamp = "1"> <Dataset name = '$dataset' interface = "default" > <Filter name = "go_parent_term" value = '$goterms'/> <Attribute name = "peptide" /> <Attribute name = "ensembl_peptide_id" /> <Attribute name = "ensembl_gene_id" /> <Attribute name = "ensembl_transcript_id" /> <Attribute name = "gene_biotype" /> <Attribute name = "external_gene_name" /> <Attribute name = "description" /> </Dataset> </Query>'
}

run_OMA()
{
    mkdir -p OMA
    cp data/parameters.drw OMA
    mkdir -p OMA/DB
    mv $1 OMA/DB/$1
    mv $2 OMA/DB/$2

    cd OMA
    export PATH=$PATH:/usr/local/OMA/bin
    OMA
}

convertsecs()
{
    ((h=${1}/3600))
    ((m=(${1}%3600)/60))
    ((s=${1}%60))
    printf "Elapsed time for total process: %02d:%02d:%02d\n" $h $m $s
}

start=$SECONDS
currPath=$(pwd)
species1=$1
species2=$2

#parse species name into ensembl dataset
dataset1=$(python get_dataset.py $species1)
dataset2=$(python get_dataset.py $species2)
goterms=$3
queryid=$4

#download appropriate ensembl sequences for species
download_seqs $dataset1 $goterms $species1\_$queryid.fa &
download_seqs $dataset2 $goterms $species2\_$queryid.fa

wait
run_OMA $species1\_$queryid.fa $species2\_$queryid.fa

#calculate ka/ks ratio for each pair of orthologous genes
cd $currPath
bash kaks_calc.sh $currPath/OMA/Output/OrthologousGroupsFasta

#export significant ka/ks pair data to txt files
mkdir -p kaksData
mkdir -p filteredGenes
python export_data.py $species1 $species2 $queryid

rm -rf OMA
rm -rf output

end=$SECONDS
elapsed=$(( end - start ))
convertsecs $elapsed