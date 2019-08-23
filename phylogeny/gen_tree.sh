#!/bin/bash
#Must have kalign and wget installed.

#Run using: ./gen_tree.sh [GO Term ID] [lifespan measurement] [optional file id]

get_data()
{

    goterms=$1
    if [ $# -lt 3 ] ; then
        id="${goterms:3}"
    else
        id=$3
    fi

    if [[ "$2" == "nla" ]] ; then
        map=data/map_nla.txt
    elif [[ "$2" == "nln1" ]] ; then
        map=data/map_nln1.txt
    else
        map=data/map_ml.txt
    fi

    #Downloads and prepares data for MSA
    download_genesets $goterms
    format_files $map

    #Performs MSA using kalign and creates phylogenetic tree with MSA using FastTree
    kalign -in sequences.fa -format fasta -out msa.afa
    ./trimal -in msa.afa -out msa.afa -gappyout
    cut -d ' ' -f 1 < msa.afa > temp.txt
    mv temp.txt msa.afa

    ./FastTree -quiet msa.afa > tree.nwk
    python tree_to_png.py $2
    
    mkdir -p trees
    mkdir -p imgs
    mv tree.nwk trees/tree\_$id.nwk
    mv tree.png imgs/tree\_$id.png
    
    rm sequences.fa
    rm msa.afa

}

download_genesets()
{

    echo "Downloading sequence data..."

    mkdir -p seqs
    mkdir -p debug
    input="../data/dataSpecies.txt"

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

}

format_files()
{
    
    map=$1
    #Create map of species name to lifespan value
    while read -r line; do 
        declare "$line" 
    done < $map

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

        numGenes=$(grep -c "^>" $file)
        
        grep -v "^>" $file | awk -v id="$species--${!species}--$numGenes" 'BEGIN { ORS=""; print ">"id"\n" } { print }' | tr -d '*' >> sequences.fa
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

get_data $1 $2 $3

end=$SECONDS
elapsed=$(( end - start ))
convert_secs $elapsed
