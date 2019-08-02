#!/bin/bash

#Prerequisites:
#Directory containing orthologue gene pairs (fasta format)
#MUSCLE unix executable, KaKs unix executable, pal2nal perl script, parseFastaIntoAxt perl script all in same directory
#EMBOSS installed

#Run using: bash kaks_calc.sh [path_to_directory]

loop_kaks()
{
    RED='\033[0;31m'
    NC='\033[0m'

    currPath=$(pwd)
    cd $1
    count=$(ls | grep -c '\.fa$')
    cd $currPath
    counter=1
    for file in $1/*; do 
        if [ -f "$file" ]; then 

            # name output file
            a=$file
            xpath=${a%/*}
            xbase=${a##*/}
            xpref=${xbase%.*}
            outFile=out_$xpref.kaks

            if [ ! -e output/$outFile ]; then 
                echo -e "${RED}Working on group" $counter "out of" $count".${NC}"
                run_kaks $file $outFile
            else
                echo "Skipping group" $counter". Already found in output directory."
            fi

            let counter++
            echo; echo;

        fi 
    done

    echo -e "${RED}FINISHED.${NC}"

}

run_kaks()
{
    # main pipeline
    backtranseq -sequence $1 -outfile dna.fa >/dev/null
    ./muscle -in $1 -out aln.fa >/dev/null
    perl pal2nal.pl aln.fa dna.fa -output fasta -nogap > nal.fa
    perl parseFastaIntoAXT.pl nal.fa >/dev/null
    ./KaKs -i nal.fa.axt -o $2
    
    # move out file to output folder
    mv $2 output/$2

    # remove temporary pipeline files
    rm dna.fa aln.fa nal.fa nal.fa.axt
}

convertsecs()
{

    ((h=${1}/3600))
    ((m=(${1}%3600)/60))
    ((s=${1}%60))
    printf "Elapsed time: %02d:%02d:%02d\n" $h $m $s

}

start=$SECONDS
mkdir -p output
if [ $# -eq 1 ]; then
    if [ -d $1 ]; then
        loop_kaks $1
        # run_kaks $1/OG1.fa out_OG1.kaks
    else
        echo "Error: Specified folder path does not exist."
    fi
else
    echo "Error: The input for the program is just the filename of the orthologue pairs of amino acid sequences."
fi

end=$SECONDS
elapsed=$(( end - start ))
convertsecs $elapsed