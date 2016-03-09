## bamProcessing

# @importFrom stringr str_locate_all str_match
# @importFrom Rsamtools asBam
#.createBamFile <- function(strings, bamFile = NULL) {
    
    ## we assume all reads have been aligned to the same reference, so use the header from the first one
#    nl.idx <- str_locate_all(strings[1], "\n")[[1]][,1]
#    samHeader <- str_sub(strings[1], start = 1, end = sort(nl.idx, decreasing = TRUE)[2]-1)
    
    ## find the entry for the read, ignore the header
#    samRecords <- str_match(strings, "\n([[:print:]\t]*)\n$")[,2]
    
#    tmpFile <- tempfile()
#    writeLines(c(samHeader, samRecords), con = tmpFile)
#    asBam(file = tmpFile, destination = bamFile)
#}