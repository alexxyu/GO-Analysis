#!/bin/bash
#To use: bash screen.sh [GO term file] [Screened species name] [Lower gene bound] [Upper gene bound]
#GO term file must be list of GO terms to screen separated by new line.
#Program keeps banlist so that results are cached and skipped if run again.

#Must have wget installed.

screen()
{
    input=$1
    species=$2
    lowerBound=$3
    upperBound=$4

    dataset=$(python get_dataset.py $species)

    while read -r term; do

        if grep -q $term data/banList.txt; then
            continue
        else

            echo "Working on GO term" $term"."
            download_genes $dataset $term $species
            numTerms=$(wc -l count/$species.fa | awk '{print $1}')
            
            if [ $numTerms -gt $upperBound ]; then
                #Remove ancestor terms since they will have more terms
                echo "Removing ancestor terms."
                Rscript filter.R $term 1
            elif [ $numTerms -lt $lowerBound ]; then
                #Remove descendant terms since they will have fewer terms
                echo "Removing descendant terms."
                Rscript filter.R $term 0
            else
                echo "${term}" >> data/filteredTerms.txt
            fi

            echo "${term}" >> data/banList.txt

        fi
    
    done < "$input"
}

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

convert_secs()
{
    ((h=${1}/3600))
    ((m=(${1}%3600)/60))
    ((s=${1}%60))
    printf "Elapsed time for total process: %02d:%02d:%02d\n" $h $m $s
}

if [ $# -lt 4 ]; then
    echo "Error -- Program takes arguments as follows: "
    echo "bash screen.sh [GO term file] [Screened species name] [Lower gene bound] [Upper gene bound]"
    exit 1
fi

start=$SECONDS

mkdir -p count
touch banList.txt
mv banList.txt data/banList.txt

screen $1 $2 $3 $4

rm -rf count

end=$SECONDS
elapsed=$(( end - start ))
convert_secs $elapsed