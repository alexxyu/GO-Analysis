#!/usr/bin/env Rscript
suppressMessages(library(GO.db))

#Outputs terms that are both in the given set of terms and the data file
outputIntersection <- function(set) {

    con <- file('data/go_terms.txt', open = 'r')
    while(TRUE) {
        line <- readLines(con, n = 1)
        if(length(line) == 0) break
        else if(line %in% set){
            write(line, file = "data/banList.txt", append = TRUE)
        } 
    }
    
}

args <- commandArgs(trailingOnly=TRUE)

#Makes sures enough arguments are supplied
if (length(args)<2) {
  stop("Two arguments must be supplied", call.=FALSE)
}

goTerm <- args[1]
func <- args[2]

#Second argument determines whether offspring or ancestor terms are filtered
#0: offspring filtered out
#1: ancestors filtered out
if (func == 0) {
    offspring = GOBPOFFSPRING[[goTerm]]
    outputIntersection(offspring)
} else {
    ancestors = GOBPANCESTOR[[goTerm]]
    outputIntersection(ancestors)
}