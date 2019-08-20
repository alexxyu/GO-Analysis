# GO Analysis Package

This package includes various tools that can be used to retrieve and analyze genes of animal species from Ensembl using Gene Ontology (GO) IDs. Tools include testing for postive selection by calculating Ka/Ks ratio for a set of gene sequences for a pair of species as well as retrieving the number of genes under a given GO ID for a given species.

## Getting Started

Download the project files and make sure prequisites are met.

### Prerequisites

* [wget](https://www.gnu.org/software/wget/) must be installed for sequence and gene data retrieval from Ensembl.
* [OMA Standalone](https://omabrowser.org/standalone/) must be downloaded and installed for all-vs-all sequence comparison for Ka/Ks ratio testing.
* [EMBOSS](http://emboss.sourceforge.net/apps/#list) must be installed.
* [R](https://www.r-project.org/) must be installed along with [Bioconductor](https://bioconductor.org/install/) and [GO.db] (https://bioconductor.org/packages/release/data/annotation/html/GO.db.html) libraries.
* Python 3.6 installed w/ numpy, pandas, goatools, scipy.

## Testing for Positive Selection

### Instructions

From project folder, cd to calcKaKs directory. Then run gen_kaks_data.sh. Parameters taken are formatted name of first species, formatted name of second species, GO term(s) separated by commas, and a unique query ID. As reference, you can find the correct formatted names of your species of interest in dataSpecies.txt in the data folder.

For example:
```
bash gen_kaks_data.sh GuineaPig Opossum GO:0004407,GO:1901727 histone_deacetylase
```

Runtime depends heavily on the number of genes under the given GO term(s) and scales dramatically as the number increases. The script will output two files, one in the directory kaksData and the other in the directory filteredGenes. Both files are named accordingly as follows: (species1)_(species2)_(query_id).txt. The example above will, for instance, yield the filename GuineaPig_Opossum_histone_deacetylase.txt.

The file in kaksData outputs raw Ka/Ks ratio data for each orthologous gene pair between the two species. The file in filteredGenes lists formatted data with genes that have tested for positive selection.

## Counting Genes for GO Correlation Analysis

### Counting Genes Under One GO Term for All Species

Run the countGO.sh script. Example:
```
bash countGO.sh GO:1902991,GO:1902992 test
```
This script will output a file in the countData directory with the filename pattern count_(query ID).txt. The output is a list of the gene counts, with each line referring to the dataSpecies.txt file in the same order.

### GO Correlation Analysis

#### Screening

Screening first with one species, rather than a whole set of animals, can dramatically lower runtime in the end if using a large set of GO terms.
```
bash screen.sh [GO term file] [Screened species name] [Lower gene bound] [Upper gene bound]
bash screen.sh data/go_terms.txt Microbat 10 120
```
This will generate an output file in data/filteredTerms.txt that can be used for counting genes in the next step.

#### Counting Genes Under Each GO Term for All Species

Run the countAllGO.sh script. The two arguments that it takes are the path to the file containing GO terms to count and the output filename.
```
bash countAllGO.sh data/filteredTerms.txt sigGOTerms.txt
```
This will generate a directory of count data. For correlation analysis, it calls GOAnalysis.py, which calculates normalized lifespan based on a linear regression model between BMR/M and lifespan and then calculates Spearman's rank correlation coefficient with p-value. GOAnalysis.py will output an entire tab-delimited text file with counts, which can be analyzed in Python, R, etc.

## Other Dependencies

* [MUSCLE](https://www.drive5.com/muscle/) - Multiple sequence alignment tool
* [PAL2NAL](http://www.bork.embl.de/pal2nal/) - Conversion tool for protein MSA to DNA sequences
* [Ka/Ks Calculator](https://code.google.com/archive/p/kaks-calculator/)
