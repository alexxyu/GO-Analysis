#!/bin/bash

get_data()
{

    goterms=$1
    id=${goterms:3}

    #Downloads and prepares data for MSA
    #download_index $goterms
    download_genesets $goterms
    format_files

    #Performs MSA using kalign and creates phylogenetic tree with MSA using Muscle
    kalign -in sequences.fa -format fasta -out msa.afa
    ./trimal -in msa.afa -out msa.afa -gappyout
    ./muscle -maketree -in msa.afa -out tree\_$id.phy

}

download_index()
{

    echo "Downloading index data..."
    echo "Complete."; echo ""

    mkdir -p seqs
    mkdir -p debug
    goterms=$1
    goterms=\"$goterms\"

    #Downloads list of human genes used as conserved list
    wget -o debug/index.txt -O seqs/index.txt --timeout=300 --tries=3 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?> <!DOCTYPE Query> <Query virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" completionStamp = "1"> <Dataset name = "hsapiens_gene_ensembl" interface = "default" > <Filter name = "go_parent_term" value = '$goterms'/> <Attribute name = "external_gene_name" /> </Dataset> </Query>' > /dev/null

}

download_genesets()
{

    echo "Downloading sequence data..."

    mkdir -p seqs
    mkdir -p debug
    input="../data/dataSpecies.txt"

    #Formats indexed gene list for data query
    #csv_string=$(paste -sd, seqs/index.txt)
    #csv_string="\"$csv_string\""
    while read -r species; do

        goterms=$1
        dataset=$(python ../get_dataset.py $species)
        download_genes $dataset $goterms

    done < "$input"

    #Count number of genes in each species datafile
    wait
    while read -r species; do
        
        #Ensures that the relevant species data file is fully downloaded
        lastLine=$(awk '/./{line=$0} END{print line}' seqs/$species.fa)
        while [[ "$lastLine" != "[success]" ]]; do

            debugLine=$(awk '/./{line=$0} END{print line}' debug/$species.txt)
            if [[ "$debugLine" == *"ERROR 500: Internal Server Error."* ]] || [[ "$debugLine" == *"unable to resolve host address"* ]] ; then
                echo "Server error thrown. Skipping current GO term."
                exit 1
            fi

            sleep 3
            lastLine=$(awk '/./{line=$0} END{print line}' seqs/$species.fa)
        done

        #Deletes "[success]" line from query
        sed -i '' -e '$ d' seqs/$species.fa

    done < "$input"
    rm -rf debug

    echo "Complete."; echo ""

}

download_genes()
{

    #Dataset name from Ensembl
    dataset=$1
    dataset=\"$dataset\"

    #GO terms must be separated by comma (no space)
    goterms=$2
    goterms=\"$goterms\"

    wget -b -o debug/$species.txt -O seqs/$species.fa --timeout=300 --tries=3 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?> <!DOCTYPE Query> <Query virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" completionStamp = "1"> <Dataset name = '$dataset' interface = "default" > <Filter name = "go_parent_term" value = '$goterms'/> <Attribute name = "peptide" /> <Attribute name = "ensembl_peptide_id" /> </Dataset> </Query>' > /dev/null
    #if [ $dataset == "\"hsapiens_gene_ensembl\"" ]; then
    #    wget -b -o debug/$species.txt -O seqs/$species.fa --timeout=300 --tries=3 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?> <!DOCTYPE Query> <Query virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" completionStamp = "1"> <Dataset name = "hsapiens_gene_ensembl" interface = "default" > <Filter name = "go_parent_term" value = '$goterms'/> <Attribute name = "peptide" /> <Attribute name = "ensembl_peptide_id" /> </Dataset> </Query>' > /dev/null
    #else
        #wget -b -o debug/$species.txt -O seqs/$species.fa --timeout=300 --tries=3 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?> <!DOCTYPE Query> <Query virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" completionStamp = "1"> <Dataset name = '$dataset' interface = "default" > <Filter name = "external_gene_name" value = '$csv_string'/> <Attribute name = "peptide" /> <Attribute name = "external_gene_name" /> </Dataset> </Query>' > /dev/null
    #fi

}

format_files()
{
    
    for file in seqs/*.fa; do
        #Removes any unavailable sequences
        grep -v "Sequence unavailable" $file > temp.fa
        mv temp.fa $file

        species="${file%.*}"
        species="${species:5}"

        #Sorts sequences based on header
        cat $file | grep "^>" | sort | while read ID ; do awk 'BEGIN{RS=">"; ORS="";} /^'${ID:1}'/{print ">" $0; exit(0);}' $file >> seqs/$species\_sorted.fa ; done
    done

    > sequences.fa
    #Concatenates sequences of same species together, then outputs these mega-sequences to master sequences file for MSA
    for file in seqs/*_sorted.fa; do
        species="${file%.*}"
        species="${species:5}"
        species="$(echo $species | cut -f1 -d"_")"
        
        grep -v "^>" $file | awk -v id="$species" 'BEGIN { ORS=""; print ">"id"\n" } { print }' | tr -d '*' >> sequences.fa
        echo "" >> sequences.fa
        echo "" >> sequences.fa
    done

    sed -i -e 's/\[success\]//g' sequences.fa
    rm sequences.fa-e
    rm -rf seqs

}

convert_secs()
{
    ((h=${1}/3600))
    ((m=(${1}%3600)/60))
    ((s=${1}%60))
    printf "Elapsed time: %02d:%02d:%02d\n" $h $m $s
}

start=$SECONDS

get_data $1

end=$SECONDS
elapsed=$(( end - start ))
convert_secs $elapsed
