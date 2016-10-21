
#' @describeIn Fast5Summary Subset object and return an object of the same class.
#' 
#' @param x Object of class Fast5Summary
#' @param i Vector defining index to subset by.
#' 
#' @docType methods
#' @export
#' 
#' @importFrom dplyr select
#' @importFrom tidyr gather
#' @importFrom XVector compact
#' @importFrom methods initialize
setMethod("[", c("Fast5Summary", "ANY"), function(x, i) {
    recordTable <- .matchRecords(x)[i,]
    
    ## there are several columns holding information on baseCalled and fastq
    ## we do some manipulation here to get a single index list.
    baseCalledIDX <- select(recordTable, baseCalledTemplate:baseCalledComplement) %>% 
        gather(key = component, value = idx) %>% 
        filter(!is.na(idx)) %>%
        select(idx)
    baseCalledIDX <- as.vector(baseCalledIDX[,1])
    
    fastqIDX <- select(recordTable, fastqTemplate:fastq2D) %>% 
        gather(key = component, value = idx) %>% 
        filter(!is.na(idx)) %>%
        select(idx)
    fastqIDX <- as.vector(fastqIDX[,1])
    
    initialize(x, 
                readInfo = x@readInfo[recordTable[,readInfo],],
                rawData = x@rawData[recordTable[,rawData],],
                baseCalled = x@baseCalled[baseCalledIDX,],
                fastq = compact(x@fastq[fastqIDX]) )
})


#' @importFrom ShortRead id
.idFromFASTQ <- function(fastq_obj) {
    as.integer(do.call(rbind, strsplit(as.character(ShortRead::id(fastq_obj)), "_"))[,1])
}

#' @importFrom ShortRead id
.readtypeFromFASTQ <- function(fastq_obj) {
    do.call(rbind, strsplit(as.character(ShortRead::id(fastq_obj)), "_"))[,2]
}

.matchRecords <- function(summaryData) {
    
    ids <- summaryData@readInfo[,id]
    ri_row <- match(ids, summaryData@readInfo[,id])
    rd_row <- match(ids, summaryData@rawData[,id])
    bct_row <- match(ids, summaryData@baseCalled[strand == "template",id])
    bcc_row <- match(ids, summaryData@baseCalled[strand == "complement",id]) + nrow(summaryData@baseCalled[strand == "template"])
    fastq_id <- .idFromFASTQ( fastq(summaryData) )
    fastq_id <- fastq_id[-(1:nrow(summaryData@baseCalled))]
    fq_row <- match(ids, fastq_id) + nrow(summaryData@baseCalled)
    
    record_table <- tibble(id = ids, 
                           readInfo = ri_row,
                           rawData = rd_row,
                           baseCalledTemplate = bct_row,
                           baseCalledComplement = bcc_row,
                           fastqTemplate = bct_row,
                           fastqComplement = bcc_row,
                           fastq2D = fq_row)
    
    return(record_table)
    
}

.get2D <- function(summaryData) {
    
    ids <- filter(baseCalled(summaryData), strand == "template", full_2D == TRUE)[,id]
    idx <- which(readInfo(summaryData)[,id] %in% ids)
    return(summaryData[idx,])
}


#' Extract template reads
#' 
#' This generic function accesses the fastq slot stored in an object derived 
#' from the Fast5Summary class, and returns only the template reads.
#' 
#' @param x Object of class \code{\linkS4class{Fast5Summary}}
#' @return A ShortReadQ object
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    fastqTemplate( s.typhi.rep2 )
#' }
setGeneric("fastqTemplate", function(x) {
    standardGeneric("fastqTemplate")
})

#' Extract complement reads
#' 
#' This generic function accesses the fastq slot stored in an object derived 
#' from the Fast5Summary class, and returns only the complement reads.
#' 
#' @param x Object of class \code{\linkS4class{Fast5Summary}}
#' @return A ShortReadQ object
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    fastqComplement( s.typhi.rep2 )
#' }
setGeneric("fastqComplement", function(x) {
    standardGeneric("fastqComplement")
})

#' Extract 2D reads
#' 
#' This generic function accesses the fastq slot stored in an object derived 
#' from the Fast5Summary class, and returns only the 2D reads.
#' 
#' @param x Object of class \code{\linkS4class{Fast5Summary}}
#' @return A ShortReadQ object
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    fastq2D( s.typhi.rep2 )
#' }
setGeneric("fastq2D", function(x) {
    standardGeneric("fastq2D")
})

#' @describeIn Fast5Summary Returns ShortReadQ object containing only template reads
#' 
#' @export
setMethod("fastqTemplate", 
          c(x = "Fast5Summary"),
          function(x) {
              idx <- which(.readtypeFromFASTQ( fastq(x) ) == 'template')
              return( fastq(x)[idx,] )
          }
)

#' @describeIn Fast5Summary Returns ShortReadQ object containing only complement reads
#' 
#' @export
setMethod("fastqComplement", 
        c(x = "Fast5Summary"),
        function(x) {
            idx <- which(.readtypeFromFASTQ( fastq(x) ) == 'complement')
            return( fastq(x)[idx,] )
        }
)

#' @describeIn Fast5Summary Returns ShortReadQ object containing only 2D reads
#' 
#' @export
setMethod("fastq2D", 
          c(x = "Fast5Summary"),
          function(x) {
              idx <- which(.readtypeFromFASTQ( fastq(x) ) == '2D')
              return( fastq(x)[idx,] )
          }
)