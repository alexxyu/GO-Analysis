#!/bin/bash
#Must have kalign and wget installed.

get_data()
{

    if [ -z "$o_flag" ] && [ ! -z "$g_flag" ] ; then
        id="${g_flag:3:7}"
    else
        id=$o_flag
    fi

    echo "Output tree filename: tree_"$id".nwk" 

    if [[ "$m_flag" == "nla" ]] ; then
        map=data/map_nla.txt
    elif [[ "$m_flag" == "nln1" ]] ; then
        map=data/map_nln1.txt
    else
        map=data/map_ml.txt
    fi

    #Downloads and prepares data for MSA
    download_genesets $g_flag
    format_files $map

    #Performs MSA using kalign and creates phylogenetic tree with MSA using FastTree
    kalign -in sequences.fa -format fasta -out msa.afa
    ./trimal -in msa.afa -out msa.afa -gappyout
    cut -d ' ' -f 1 < msa.afa > temp.txt
    mv temp.txt msa.afa

    ./FastTree -quiet msa.afa > tree.nwk
    python tree_to_png.py $d_flag $m_flag
    
    mkdir -p trees
    mkdir -p imgs
    mv tree.nwk trees/tree\_$id.nwk
    mv tree.png imgs/tree\_$id.png
    
    rm sequences.fa
    rm msa.afa

}

download_genesets()
{

    mkdir -p seqs
    mkdir -p debug
    input="../data/dataSpecies.txt"

    if [ ! -z "$g_flag" ]; then
        echo "Downloading GO term sequence data..."
        goterms=$1
        goterms=\"$goterms\"
    else
        echo "Downloading specified gene sequence data..."
        geneList=$(tr -s '\n' ',' < $l_flag | sed 's/.$//')
        geneList=\"$geneList\"
    fi

    while read -r species; do

        dataset=$(python ../get_dataset.py $species)

        if [ ! -z "$g_flag" ]; then
            download_genes $dataset $goterms
        else
            download_genes_by_list $dataset $geneList
        fi

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

download_genes_by_list()
{

    #Dataset name from Ensembl
    dataset=$1
    dataset=\"$dataset\"

    geneList=$2

    wget -b -o debug/$species.txt -O seqs/$species.fa --timeout=300 --tries=3 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?> <!DOCTYPE Query> <Query virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" completionStamp = "1"> <Dataset name = '$dataset' interface = "default" > <Filter name = "external_gene_name" value = '$geneList'/> <Attribute name = "peptide" /> <Attribute name = "external_gene_name" /> </Dataset> </Query>' > /dev/null

}

download_genes()
{

    #Dataset name from Ensembl
    dataset=$1
    dataset=\"$dataset\"

    goterms=$2

    wget -b -o debug/$species.txt -O seqs/$species.fa --timeout=300 --tries=3 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?> <!DOCTYPE Query> <Query virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" completionStamp = "1"> <Dataset name = '$dataset' interface = "default" > <Filter name = "go_parent_term" value = '$goterms'/> <Attribute name = "peptide" /> <Attribute name = "external_gene_name" /> </Dataset> </Query>' > /dev/null

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

        #Deletes duplicate genes
        if [[ "$r_flag" == "true" ]]; then
            awk '/^>/{f=!d[$1];d[$1]=1}f' seqs/$species\_sorted.fa > temp.fa
            mv temp.fa seqs/$species\_sorted.fa
        fi
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

print_usage() {
    echo "Usage: -d to disable color view in phylogeny tree"
    echo "       -g to specify GO terms to use (overrides -l)"
    echo "       -l to specify file from which to use gene list"
    echo "       -m to specify lifespan measurement (default: maximum lifespan)"
    echo "       -o to specify output filename id"
    echo "       -r to remove duplicate sequences"
    echo ""
}

validate_inputs() {

    if [ ! -z "$l_flag" ] && [ ! -f "$l_flag" ]; then
        echo "Specified gene list file does not exist."
        exit 1
    fi

    if [ -z "$g_flag" ] && [ -z "$l_flag" ]; then
        echo "Have to specify either GO term list or gene list filename."
        exit 1
    fi

}

start=$SECONDS

d_flag='true'   # disable color view in phylogeny tree
g_flag=''       # specify GO terms to use
l_flag=''       # specify gene list file to use
m_flag='ml'     # specify lifespan measurement
o_flag=''       # specify output filename id
r_flag='false'  # remove duplicate genes

while getopts 'dg:l:m:ro:' flag; do
  case "${flag}" in
    d) d_flag='false' ;;
    g) g_flag="${OPTARG}" ;;
    l) l_flag="${OPTARG}" ;;
    m) m_flag="${OPTARG}" ;;
    r) r_flag='true' ;;
    o) o_flag="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

validate_inputs
get_data

end=$SECONDS
elapsed=$(( end - start ))
convert_secs $elapsed
