#' Read summary data from fast5 files.
#' 
#' Reads one or more fast5 files and collects summary information about them.  
#' 
#' Currently this function assumes all files passed to it come from the same 
#' sequencing run.  It makes no effort to check for alternative file names or 
#' the like.  If files from multiple runs are passed to it they will be 
#' collated together and any analysis performed on them will reprsent the 
#' mixture of both experiments.
#' 
#' @param files Character vector of fast5 files to be read.
#' @return Object of class \linkS4class{Fast5Summary}
#' @examples \dontrun{
#' fast5files <- list.file('/foo/bar/', pattern = '.fast5$')
#' summaryData <- readFast5Summary(fast5files)
#' }
#' @export
#' @importFrom dplyr filter
#' @importFrom ShortRead append
readFast5Summary <- function(files) {
    
    message("Reading Channel Data")
    readInfo <- lapply(files, .getReadChannelMux)
    readInfo <- rbindlist(readInfo)
    readInfo <- data.table(id = 1:nrow(readInfo), file = basename(files), readInfo)
    
    message("Reading Raw Data")
    rawData <- lapply(files, .getSummaryRaw)
    rawData <- data.table(id = readInfo[,id], rbindlist(rawData))
    
    message("Reading Template Data")
    template <- lapply(files, .getSummaryBaseCalled, strand = "template")
    template <- data.table(id = readInfo[,id], rbindlist(template))
    template <- filter(template, !(is.na(num_events)))
    
    message("Reading Complement Data")
    complement <- lapply(files, .getSummaryBaseCalled, strand = "complement")
    complement <- data.table(id = readInfo[,id], rbindlist(complement))
    complement <- filter(complement, !(is.na(num_events)))
    
    message("Reading Template FASTQ")
    ## get the fastq for those that have it
    fq_t <- sapply(files[ template[,id] ], .getFastqString, strand = "template")
    fq_t <- .processFastqVec(fq_t, readIDs = template[, id], appendID = "_template")  
    fastq_template <- fq_t$fastq
    ## if there are any invalid entries we need to remove them
    if(length(fq_t$invalid)) {
        template <- template[-fq_t$invalid,]
    }
    
    message("Reading Complement FASTQ")
    fq_c <- sapply(files[ complement[,id] ], .getFastqString, strand = "complement")
    fq_c <- .processFastqVec(fq_c, readIDs = complement[, id], appendID = "_complement")  
    fastq_complement <- fq_c$fastq
    ## if there are any invalid entries we need to remove them
    if(length(fq_c$invalid)) {
        complement <- complement[-fq_c$invalid,]
    }
    
    message("Reading 2D FASTQ")
    ## we haven't read anything about 2D reads yet, so we need to identify which
    ## files have them.  Then we'll read only those
    idx2D <- which(sapply(files, .groupExistsString, group = paste0("/Analyses/Basecall_2D_000/BaseCalled_2D")))
    fq_2D <- sapply(files[ idx2D ], .getFastqString, strand = "2D")
    fq_2D <- .processFastqVec(fq_2D, readIDs = readInfo[idx2D, id], appendID = "_2D")  
    fastq_2D <- fq_2D$fastq
    ## if there are any invalid entries we need to remove them
    if(length(fq_2D$invalid)) {
        twoD <- twoD[-fq_c$invalid,]
    }
    
    ## We update the individual strands to indicate if they are part of a full 2D read
    template <- data.table(template, full_2D = template[,id] %in% idx2D)
    complement <- data.table(complement, full_2D = complement[,id] %in% idx2D)
    
    ## combine the template, complement and 2D data
    baseCalled <- rbind(template, complement)
    fastq <- append(fastq_template, fastq_complement)
    fastq <- append(fastq, fastq_2D)
    
    message("Done")
    obj <- new("Fast5Summary", readInfo = readInfo, rawData = rawData, baseCalled = baseCalled, fastq = fastq)
    
    return(obj)
}