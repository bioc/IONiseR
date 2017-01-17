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
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr filter
#' @importFrom ShortRead append
#' @importFrom data.table rbindlist
readFast5Summary <- function(files) {
    
    ## some files can't be opened, so we filter them here
    message("Checking file validity")
    fileStatus <- sapply(files, .checkOpening, USE.NAMES = FALSE)
    files <- files[ which(fileStatus) ]
    
    ## use a subset of files to try and guess the version
    versions <- sapply(sample(files, size = min(length(files), 15)), .fast5version, USE.NAMES = FALSE)
    if(all(versions) - mean(versions) == 0) {
        warning("Not all files appear to be processed with the same version of MinKNOW.\nThis may cause problems later")
    }
    
    message("Reading Channel Data")
    readInfo <- lapply(files, .getReadChannelMux)
    readInfo <- rbindlist(readInfo)
    readInfo <- as_tibble(cbind(id = 1:nrow(readInfo), file = basename(files), readInfo))
    ## if the files don't meet our expected structure, we should catch it here
    if(all(is.na(readInfo[['read']]))) {
        stop("No files matched the expected fast5 file structure.\n  Have they been processed and basecalled with MinKNOW?")
    } else if (length(idx <- which(is.na(readInfo[['read']])))) {
        message(length(idx), " file(s) did not match the expected structure and have been removed.")
        readInfo <- readInfo[ -idx, ]
        readInfo[['id']] <- 1:nrow(readInfo)
        files <- files[ -idx ]
    }
    
    message("Reading Event Data")
    eventData <- lapply(files, .getSummaryEvents)
    eventData <- as_tibble(cbind(readInfo['id'], rbindlist(eventData)))
    ## we convert timing data into seconds. 
    ## To do this we find the sampling rate stored in one file
    samplingRate <- .getSamplingRate(files[1])
    eventData <- mutate(eventData, 
                      start_time = start_time / samplingRate,
                      duration = duration / samplingRate)
    
    message("Reading Template Data")
    template <- lapply(files, .getSummaryBaseCalled, strand = "template")
    template <- as_tibble(cbind(readInfo['id'], rbindlist(template)))
    template <- filter(template, !(is.na(num_events)))
    
    message("Reading Complement Data")
    complement <- lapply(files, .getSummaryBaseCalled, strand = "complement")
    complement <- as_tibble(cbind(readInfo['id'], rbindlist(complement)))
    complement <- filter(complement, !(is.na(num_events)))
    
    message("Reading Template FASTQ")
    ## get the fastq for those that have it
    fq_t <- sapply(files[ template[['id']] ], .getFastqString, strand = "template")
    fq_t <- .processFastqVec(fq_t, readIDs = template[['id']], appendID = "_template")  
    fastq_template <- fq_t$fastq
    ## if there are any invalid entries we need to remove them
    if(length(fq_t$invalid)) {
        template <- template[-fq_t$invalid,]
    }
    
    message("Reading Complement FASTQ")
    fq_c <- sapply(files[ complement[['id']] ], .getFastqString, strand = "complement")
    fq_c <- .processFastqVec(fq_c, readIDs = complement[['id']], appendID = "_complement")  
    fastq_complement <- fq_c$fastq
    ## if there are any invalid entries we need to remove them
    if(length(fq_c$invalid)) {
        complement <- complement[-fq_c$invalid,]
    }
    
    message("Reading 2D FASTQ")
    ## we haven't read anything about 2D reads yet, so we need to identify which
    ## files have them.  Then we'll read only those
    idx2D <- which(sapply(files, .groupExistsString, group = paste0("/Analyses/Basecall_2D_000/BaseCalled_2D")))
    if(length(idx2D)) {
        fq_2D <- sapply(files[ idx2D ], .getFastqString, strand = "2D")
        fq_2D <- .processFastqVec(fq_2D, readIDs = readInfo[['id']][idx2D], appendID = "_2D")  
        fastq_2D <- fq_2D$fastq
    }

    ## We update the individual strands to indicate if they are part of a full 2D read
    template <- as_tibble(cbind(template, full_2D = template[['id']] %in% idx2D))
    complement <- as_tibble(cbind(complement, full_2D = complement[['id']] %in% idx2D))
    
    ## combine the template, complement and 2D data
    baseCalled <- rbind(template, complement)
    fastq <- ShortRead::append(fastq_template, fastq_complement)
    if(length(idx2D)) {
        fastq <- ShortRead::append(fastq, fastq_2D)
    }
    
    message("Done")
    obj <- new("Fast5Summary", 
               readInfo = readInfo, 
               rawData = tibble(), 
               eventData = eventData, 
               baseCalled = baseCalled, 
               fastq = fastq,
               versions = list('IONiseR' = strsplit(as.character(packageVersion("IONiseR")),".",fixed=T)[[1]],
                               'MinKNOW' = max(versions))
               )
    
    return(obj)
}



readFast5Summary_update <- function(files) {
    
    ## some files can't be opened, so we filter them here
    message("Checking file validity")
    fileStatus <- sapply(files, .checkOpening, USE.NAMES = FALSE)
    files <- files[ which(fileStatus) ]
    
    ## use a subset of files to try and guess the version
    #versions <- sapply(sample(files, size = min(length(files), 15)), .fast5version, USE.NAMES = FALSE)
    #if(all(versions) - mean(versions) == 0) {
    #    warning("Not all files appear to be processed with the same version of MinKNOW.\nThis may cause problems later")
    #}
    
    status <- .fast5status(files = sample(files, size = min(length(files), 15)))
    
    message("Reading Channel Data")
    if(status$read_in_name) {
        readNums <- as.integer(str_match(string = files, pattern = "read([0-9]+)")[,2])
    } else { 
        readNums <- NULL 
    }
    readInfo <- do.call("rbind", mapply(.getReadChannelMux2, files, readNums, dontCheck = TRUE,
                                        USE.NAMES = FALSE, SIMPLIFY = FALSE))
    readInfo <- as_tibble(cbind(id = 1:nrow(readInfo), file = basename(files), readInfo))

    message("Reading Raw Data")
    rawData <- do.call("rbind", mapply(.getRawSummary, files, readNums, dontCheck = TRUE, 
                                       USE.NAMES = FALSE, SIMPLIFY = FALSE))

    message("Reading Event Data")
    eventData <- do.call("rbind", mapply(.getEventsSummary,  files, readNums, dontCheck = TRUE, 
                                         USE.NAMES = FALSE, SIMPLIFY = FALSE))
    rawEventData <- as_tibble(cbind(id = readInfo[['id']], rawData, eventData))
    
    
    ## we convert timing data into seconds. 
    ## To do this we find the sampling rate stored in one file
    samplingRate <- IONiseR:::.getSamplingRate(files[1])
    rawEventData <- mutate(rawEventData, 
                      start_time = start_time / samplingRate,
                      duration = duration / samplingRate)

    # message("Reading Template Data")
    # template <- lapply(files, .getSummaryBaseCalled, strand = "template")
    # template <- as_tibble(cbind(readInfo['id'], rbindlist(template)))
    # template <- filter(template, !(is.na(num_events)))
    # 
    # message("Reading Complement Data")
    # complement <- lapply(files, .getSummaryBaseCalled, strand = "complement")
    # complement <- as_tibble(cbind(readInfo['id'], rbindlist(complement)))
    # complement <- filter(complement, !(is.na(num_events)))
    # 
    # message("Reading Template FASTQ")
    # ## get the fastq for those that have it
    # fq_t <- sapply(files[ template[['id']] ], .getFastqString, strand = "template")
    # fq_t <- .processFastqVec(fq_t, readIDs = template[['id']], appendID = "_template")  
    # fastq_template <- fq_t$fastq
    # ## if there are any invalid entries we need to remove them
    # if(length(fq_t$invalid)) {
    #     template <- template[-fq_t$invalid,]
    # }
    # 
    # message("Reading Complement FASTQ")
    # fq_c <- sapply(files[ complement[['id']] ], .getFastqString, strand = "complement")
    # fq_c <- .processFastqVec(fq_c, readIDs = complement[['id']], appendID = "_complement")  
    # fastq_complement <- fq_c$fastq
    # ## if there are any invalid entries we need to remove them
    # if(length(fq_c$invalid)) {
    #     complement <- complement[-fq_c$invalid,]
    # }
    # 
    # message("Reading 2D FASTQ")
    # ## we haven't read anything about 2D reads yet, so we need to identify which
    # ## files have them.  Then we'll read only those
    # idx2D <- which(sapply(files, .groupExistsString, group = paste0("/Analyses/Basecall_2D_000/BaseCalled_2D")))
    # if(length(idx2D)) {
    #     fq_2D <- sapply(files[ idx2D ], .getFastqString, strand = "2D")
    #     fq_2D <- .processFastqVec(fq_2D, readIDs = readInfo[['id']][idx2D], appendID = "_2D")  
    #     fastq_2D <- fq_2D$fastq
    # }
    # 
    # ## We update the individual strands to indicate if they are part of a full 2D read
    # template <- as_tibble(cbind(template, full_2D = template[['id']] %in% idx2D))
    # complement <- as_tibble(cbind(complement, full_2D = complement[['id']] %in% idx2D))
    # 
    # ## combine the template, complement and 2D data
    # baseCalled <- rbind(template, complement)
    # fastq <- ShortRead::append(fastq_template, fastq_complement)
    # if(length(idx2D)) {
    #     fastq <- ShortRead::append(fastq, fastq_2D)
    # }
    
    message("Done")
    obj <- new("Fast5Summary", 
               readInfo = readInfo, 
               rawData = tibble(), 
               eventData = eventData, 
               baseCalled = baseCalled, 
               fastq = fastq,
               versions = list('IONiseR' = strsplit(as.character(packageVersion("IONiseR")),".",fixed=T)[[1]],
                               'MinKNOW' = max(versions))
    )
    
    return(obj)
}
