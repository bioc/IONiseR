.fast5status <- function(files, warn = FALSE) {
    
    ## is the read number present in the file name?
    readInName <- stringr::str_detect(string = files, pattern = "read([0-9]+)")
    
    ## is /Raw/Reads present?
    rawReads <- sapply(files, IONiseR:::.groupExistsString, group = "/Raw/Reads", USE.NAMES = FALSE)
    
    ## is /Analysis/EventDetection_000 present?
    eventDetection <- sapply(files, IONiseR:::.groupExistsString, group = "/Analyses/EventDetection_000/Reads", USE.NAMES = FALSE)
    
    ## how many times has event detection been run?
    eventDetectionNum <- sapply(files, IONiseR:::.findAnalysisNumber, grepString = "EventDetection_([0-9]+)", USE.NAMES = FALSE)
    if(length(unique(eventDetectionNum)) != 1) {
        if(warn) {
            warning("Inconsistent event detection runs.  Defaulting to the earliest", call. = FALSE)
        }
        eventNum <- "000"
    } else {
        eventNum <- unique(eventDetectionNum)
    }
    
    ## is /Analysis/Basecall_1D_000 present?
    basecall_1d <- all(sapply(files, IONiseR:::.groupExistsString, group = "/Analyses/Basecall_1D_000", USE.NAMES = FALSE))
    d <- "1D"
    
    ## is /Analysis/Basecall_2D_000 present?
    basecall_2d <- all(sapply(files, IONiseR:::.groupExistsString, group = "/Analyses/Basecall_2D_000", USE.NAMES = FALSE))
    if(basecall_2d) {
        d <- "2D"
    }
    
    ## is /Analysis/Basecall_2D_000/BaseCalled_2D present?
    analysis_2d <- all(sapply(files, IONiseR:::.groupExistsString, group = "/Analyses/Basecall_2D_000/BaseCalled_2D", USE.NAMES = FALSE))
    
    baseCallingNum <- sapply(files, IONiseR:::.findAnalysisNumber, grepString = "Basecall_[12]D_([0-9]+)", USE.NAMES = FALSE)
    if(length(unique(eventDetectionNum)) != 1) {
        if(warn) {
            warning("Inconsistent base calling runs.  Defaulting to the earliest", call. = FALSE)
        }
        basecallNum <- "000"
    } else {
        basecallNum <- unique(eventDetectionNum)
    }
    
    return(list(read_in_name = all(readInName),
                raw_reads = all(rawReads),
                event_detection = all(eventDetection),
                event_num = eventNum,
                basecall_1d = basecall_1d,
                basecall_2d = basecall_2d,
                analysis_2d = analysis_2d,
                basecall_num = basecallNum,
                d = d))
    
}
