
#' An S4 class for summarised data from a MinION sequencing run
#'
#' @slot readInfo Object of class tibble.  Contains five columns: 
#' \itemize{
#'   \item id - an integer key that allows use to match entries in the separate 
#'   slots of this object.
#'   \item file - Basename of the fast5 file the data was read from.
#'   \item read - Read number from channel.
#'   \item channel - channel.
#'   \item mux - Specific pore that was used within the four that are assigned 
#'   to a single channel. Should be in the range 1-4, but if this isn't 
#'   available it will be 0.
#' }
#' @slot eventData Object of class tibble.  Holds summary of events data 
#' prior to base calling. Contains five columns: 
#' \itemize{
#'   \item id - an integer key that allows use to match entries in the 
#'   separate slots of this object.
#'   \item start_time - time in seconds after the run started that this 
#'   reading began. 
#'   \item duration - time in seconds the reading lasted. 
#'   \item num_events - the number of events that were recorded as part of 
#'   this reading. 
#'   \item median_signal - median of the recorded signals for this set of events.
#' }
#' @slot baseCalled Object of class tibble.  For the most part contains 
#' similar data to the @@eventData slot, the base called data is derived from it.
#' \itemize{
#'   \item id - an integer key that allows use to match entries in the 
#'   separate slots of this object.
#'   \item start_time - time in seconds after the run started that this 
#'   reading began. 
#'   \item duration - time in seconds the reading lasted. 
#'   \item num_events - the number of events that were recorded as part of 
#'   this reading. 
#'   \item strand - can be either 'template' or 'complement'
#'   \item full_2D - boolean value specifying whether the read forms part of a 
#'   2D pair.  If TRUE the FASTQ data for the template, complement and 2D read 
#'   will be available in the @@fastq slot.
#' }
#' @slot fastq Object of class ShortReadQ.  This slot contains all reads 
#' (template, complement and 2D).  The read names take the form NUM_STRAND, 
#' where NUM matches with the id column in the other slots and STRAND indicates 
#' whether the read is template, complement or 2D.
#' @return An object of class Fast5Summary
#' @name Fast5Summary-class
#' @exportClass Fast5Summary
setClass("Fast5Summary",
         slots = list(readInfo = "tbl_df",
                      rawData = "tbl_df",
                      eventData = "tbl_df",
                      baseCalled = "tbl_df",
                      fastq = "ShortReadQ"))



setMethod(show, "Fast5Summary", function(object) {
    cat("Object of class: Fast5Summary\nContains information from:\n")
    cat(" ", nrow(readInfo(object)), "fast5 files\n")
    cat("  |-", nrow(filter(baseCalled(object), strand == "template")), "template strands\n")
    cat("  |-", nrow(filter(baseCalled(object), strand == "complement")), "complement strands\n")
    cat("  |-", nrow(filter(baseCalled(object), strand == "template", 
                            full_2D == TRUE)), "full 2D reads\n")
    ## if we have pass/fail information, print that too
    ## we also check that the pass column exists, hopefully we can remove this at some point
    if("pass" %in% names(readInfo(object))) {
        pf <- nrow(filter(readInfo(object), !is.na(pass)))
        if(pf) {
            cat("  |-", nrow(filter(readInfo(object), pass==TRUE)), "pass reads")
        }
    }
})

#' @describeIn Fast5Summary Returns the number of files read during creation of the object
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    length( s.typhi.rep2 )
#' }
#' @export
setMethod(length, "Fast5Summary", function(x) {
    nrow( readInfo(x) )
})



