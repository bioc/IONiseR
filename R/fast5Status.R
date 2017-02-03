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

#' @importFrom stringr str_detect
#' @importFrom dplyr transmute slice n select
.strandExistence <- function(ls, strand = "BaseCalled_template") {
    
    loc <- filter(ls, name == strand) %>% 
            select(group, name) %>% 
            transmute(path = paste(group, name, sep = "/")) %>% 
            slice(n())
    loc <- ifelse(stringr::str_detect(loc, strand), loc[[1]], "")
    return(loc)
}

.chooseStrand <- function(paths) {
    matches <- str_match(string = paths, pattern = "(^.*)([12]D_)([0-9]+)(.*$)")
    matches <- matches[which(!is.na(matches[,1])),,drop=FALSE]
    if(!nrow(matches)) {
        return("")
    }
    
    ## if the we have 1D and 2D analysis workflows, throw an error.  
    ## These are really different!
    if(length(unique(matches[,3])) != 1) {
        stop("Inconsistent analysis workflows detected.  ",
             "Were these files analyses with different versions of MinKNOW?")
    }
    
    if(length(unique(matches[,4])) != 1) {
        ## if there's a mix of anaylsis numbers, warn we're picking the lowest in the group
        warning("Inconsistent analysis runs detected.  ",
                "Defaulting to the earliest", call. = FALSE)
    } 
    return( paste0(matches[1,2], matches[1,3], sort(unique(matches[,4]))[1], matches[1,5]) )
}

#' @importFrom stringr str_detect
.fast5status_2 <- function(files, warn = FALSE) {
    
    ## is the read number present in the file name?
    readInName <- stringr::str_detect(string = files, pattern = "read([0-9]+)")
    
    lsList <- lapply(files, h5ls, recursive = 3, datasetinfo = FALSE)
    
    ## is /Raw/Reads present?
    rawReads <- sapply(lsList, FUN = function(x) {
        select(x, group) %>% 
            str_detect(pattern = "/Raw/Reads")
    })
    
    ## is /Analysis/EventDetection_000 present?
    eventDetection <- sapply(lsList, FUN = function(x) {
        select(x, group) %>% 
            str_detect(pattern = "/EventDetection")
    })
    
    template_loc <- sapply(lsList, .strandExistence, strand = "BaseCalled_template")
    template_loc <- .chooseStrand(template_loc)
    
    complement_loc <- sapply(lsList, .strandExistence, strand = "BaseCalled_complement")
    complement_loc <- .chooseStrand(complement_loc)

    return(list(read_in_name = all(readInName),
                raw_reads = all(rawReads),
                event_detection = all(eventDetection),
                template_loc = template_loc,
                complement_loc = complement_loc)) 
    
}

