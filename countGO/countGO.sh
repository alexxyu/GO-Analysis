#!/bin/bash
#To use: ./countGO.sh [GO Term(s)] [unique query id]
#For multiple GO terms, enter them separated by only a comma
#Example: ./countGO.sh GO:1902991,GO:1902992 test

#Must have wget installed.

count_GO_Terms()
{
    queryid=$2

    #Download all species data in background first
    input="../data/dataSpecies.txt"
    numSpecies=$(wc -l $input | awk '{print $1}')

    while read -r species; do

        goterms=$1
        dataset=$(python ../get_dataset.py $species)
        download_genes $dataset $goterms $species

    done < "$input"

    #Count number of genes in each species datafile
    wait
    while read -r species; do
        
        #Ensures that the relevant species data file is fully downloaded
        lastLine=$(awk '/./{line=$0} END{print line}' count/$species.fa)
        while [[ "$lastLine" != "[success]" ]]; do

            debugLine=$(awk '/./{line=$0} END{print line}' debug/$species.txt)
            if [[ "$debugLine" == *"ERROR 500: Internal Server Error."* ]] || [[ "$debugLine" == *"unable to resolve host address"* ]] ; then
                echo "Server error thrown. Skipping current GO term."
                rm countData/count\_$queryid.txt
                exit 1
            fi

            sleep 3
            lastLine=$(awk '/./{line=$0} END{print line}' count/$species.fa)
        done

        numTerms=$(wc -l count/$species.fa | awk '{print $1}')
        numTerms=$(( numTerms - 1 ))
        echo "${numTerms}" >> countData/count\_$queryid.txt

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

    wget -b -o debug/$species.txt -O $species.fa --timeout=300 --tries=3 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?> <!DOCTYPE Query> <Query virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" completionStamp = "1"> <Dataset name = '$dataset' interface = "default" > <Filter name = "go_parent_term" value = '$goterms'/> <Attribute name = "ensembl_gene_id" /> </Dataset> </Query>' > /dev/null
    mv $species.fa count/$species.fa
}

validate_file()
{
    #First checks if file already exists. If so, then checks if it contains correct number of lines.
    #If not, deletes and re-calculates. If so, skips.
    file="countData/count_$1.txt"
    if [[ -e $file ]]; then

        speciesNum=$(wc -l data/dataSpecies.txt | awk '{print $1}')
        lineCount=$(wc -l $file | awk '{print $1}')
        if [[ $lineCount -ne $speciesNum ]]; then
            rm $file
            echo $file "does not contain the right number of terms. Deleted and will re-calculate."
            printf ""
        else
            exit 0
        fi

    fi

}

convert_secs()
{
    ((h=${1}/3600))
    ((m=(${1}%3600)/60))
    ((s=${1}%60))
    printf "Elapsed time: %02d:%02d:%02d\n" $h $m $s
}

if [ $# -lt 2 ]; then
    echo "Error -- Program takes arguments as follows: "
    echo "bash countGO.sh [GO Term(s)] [unique query id]"
    exit 1
fi

start=$SECONDS

mkdir -p count
mkdir -p countData
mkdir -p debug

validate_file $2
count_GO_Terms $1 $2

rm -rf count
rm -rf debug

end=$SECONDS
elapsed=$(( end - start ))
convert_secs $elapsed
