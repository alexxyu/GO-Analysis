#!/bin/bash

#Run using: bash extract_data.sh [Dataset ID] [GO terms] [outputFilename]

#Dataset name from Ensembl
dataset=$1
dataset=\"$dataset\"

#GO terms must be separated by comma (no space)
goterms=$2
goterms=\"$goterms\"

filename=$3
wget -O $3 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?> <!DOCTYPE Query> <Query virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" completionStamp = "1"> <Dataset name = '$dataset' interface = "default" > <Filter name = "go_parent_term" value = '$goterms'/> <Attribute name = "peptide" /> <Attribute name = "ensembl_peptide_id" /> <Attribute name = "ensembl_gene_id" /> <Attribute name = "ensembl_transcript_id" /> <Attribute name = "gene_biotype" /> <Attribute name = "external_gene_name" /> <Attribute name = "description" /> </Dataset> </Query>'