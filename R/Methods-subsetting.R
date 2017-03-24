
#' @describeIn Fast5Summary Subset object and return an object of the same class.
#' 
#' @param x Object of class Fast5Summary
#' @param i Vector defining index to subset by.
#' 
#' @docType methods
#' @export
#' 
#' @importFrom dplyr select slice
#' @importFrom tidyr gather
#' @importFrom XVector compact
#' @importFrom methods initialize
setMethod("[", c("Fast5Summary", "ANY"), function(x, i) {
    recordTable <- .matchRecords(x)[i,]
    
    ## there are several columns holding information on baseCalled and fastq
    ## we do some manipulation here to get a single index list.
    baseCalledIDX <- select(recordTable, baseCalledTemplate:baseCalledComplement) %>% 
        gather(key = component, value = idx) %>% 
        filter(!is.na(idx))
    baseCalledIDX <- baseCalledIDX[['idx']]
    
    fastqIDX <- select(recordTable, fastqTemplate:fastq2D) %>% 
        gather(key = component, value = idx) %>% 
        filter(!is.na(idx))
    fastqIDX <- fastqIDX[['idx']]
    
    if(.hasSlot(x, 'versions')) {
        versions <- x@versions
    } else {
        versions <- list()
    }
    
    initialize(x,
               readInfo = slice(readInfo(x), recordTable[['readInfo']]),
               eventData = slice(eventData(x), recordTable[['eventData']]),
               baseCalled = slice(baseCalled(x), baseCalledIDX),
               fastq = XVector::compact(fastq(x)[fastqIDX]),
               versions = versions)
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
    
    ids <- readInfo(summaryData)[['id']]
    ri_row <- match(ids, readInfo(summaryData)[['id']])
    ed_row <- match(ids, eventData(summaryData)[['id']])
    bct_row <- match(ids, filter(baseCalled(summaryData), strand == "template")[['id']])
    bcc_row <- match(ids, filter(baseCalled(summaryData), strand == "complement")[['id']]) + 
        nrow(filter(baseCalled(summaryData), strand == "template"))
    
    fastq_id <- .idFromFASTQ( fastq(summaryData) )
    fastq_id <- fastq_id[-(1:nrow(summaryData@baseCalled))]
    fq_row <- match(ids, fastq_id) + nrow(summaryData@baseCalled)
    
    record_table <- tibble(id = ids, 
                           readInfo = ri_row,
                           eventData = ed_row,
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