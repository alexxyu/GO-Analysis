library(GO.db, quietly = TRUE, warn.conflicts = FALSE)

outputIntersection <- function(set) {

    con <- file('data/go_terms.txt', open = 'r')
    while(TRUE) {
        line <- readLines(con, n = 1)
        if(length(line) == 0) break
        else if(line %in% set){
            write(line, file = "banList.txt", append = TRUE)
        } 
    }
    
}

args <- commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("At least two arguments must be supplied", call.=FALSE)
}

goTerm <- args[1]
func <- args[2]

if (func == 0) {
    offspring = GOBPOFFSPRING[[goTerm]]
    outputIntersection(offspring)
} else {
    ancestors = GOBPANCESTOR[[goTerm]]
    outputIntersection(ancestors)
}