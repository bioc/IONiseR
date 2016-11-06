#' Extract readInfo slot
#' 
#' This generic function accesses the readInfo slot stored in an object 
#' derived from the Fast5Summary class.
#' 
#' @param x Object of class \code{\linkS4class{Fast5Summary}}
#' @return A data.frame with 5 columns
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    readInfo( s.typhi.rep2 )
#' }
setGeneric("readInfo", function(x) {
    standardGeneric("readInfo")
})

#' @describeIn Fast5Summary Returns readInfo data.frame
#' 
#' @include classes.R
#' @export
setMethod("readInfo", 
          c(x = "Fast5Summary"),
          function(x) {
              x@readInfo
          }
)


#' Extract eventData slot
#' 
#' This generic function accesses the eventData slot stored in an object derived 
#' from the Fast5Summary class.
#' 
#' @param x Object of class \code{\linkS4class{Fast5Summary}}
#' @return A data.frame with 5 columns
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    eventData( s.typhi.rep2 )
#' }
setGeneric("eventData", function(x) {
    standardGeneric("eventData")
})

#' @describeIn Fast5Summary Returns eventData data.frame
#' 
#' @include classes.R
#' @export
setMethod("eventData", 
          c(x = "Fast5Summary"),
          function(x) {
              x@eventData
          }
)

#' Extract baseCalled slot
#' 
#' This generic function accesses the baseCalled slot stored in an object 
#' derived from the Fast5Summary class.
#' 
#' @param x Object of class \code{\linkS4class{Fast5Summary}}
#' @return A data.frame with 6 columns
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    baseCalled( s.typhi.rep2 )
#' }
setGeneric("baseCalled", function(x) {
    standardGeneric("baseCalled")
})


#' @describeIn Fast5Summary Returns baseCalled data.frame
#' 
#' @include classes.R
#' @export
setMethod("baseCalled", 
          c(x = "Fast5Summary"),
          function(x) {
              x@baseCalled
          }
)

#' Extract fastq slot
#' 
#' This generic function accesses the fastq slot stored in an object 
#' derived from the Fast5Summary class.
#' 
#' @param x Object of class \code{\linkS4class{Fast5Summary}}
#' @return A ShortReadQ object
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    fastq( s.typhi.rep2 )
#' }
setGeneric("fastq", function(x) {
    standardGeneric("fastq")
})

#' @describeIn Fast5Summary Returns ShortReadQ object stored in fastq slot.
#' 
#' @include classes.R
#' @export
setMethod("fastq", 
          c(x = "Fast5Summary"),
          function(x) {
              x@fastq
          }
)