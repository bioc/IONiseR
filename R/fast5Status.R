#' @importFrom stringr str_detect
#' @importFrom dplyr transmute slice select n
#' @noRd 
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
    matches <- matches[which(!is.na(matches[,1])),,drop = FALSE]

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
        ## if there's a mix of anaylsis numbers, warn we're picking the lowest
        warning("Inconsistent number of analysis runs detected between files.",
                " Defaulting to the earliest", call. = FALSE)
    } 
    return( paste0(matches[1,2], matches[1,3], 
                   sort(unique(matches[,4]))[1], matches[1,5]) )
}

#' Determines the 'processing status' of the files that have been supplied.
#' 
#' Fast5 files can be obtained at several points in the standard processing
#' flow.  Exactly which point was reached determines the type of data that is
#' present in the files e.g. raw signal, events detected, bases called etc.
#' This function looks at the hdf5 structure of the files to try and determine
#' what to expect.  Knowing about the structure can improve the performance of
#' data extraction code, since we don't necessarily have to check for 
#' existance every time.  This also trys to determine whether the '1D' or '2D'
#' workflow was run, since this also influences the path to (and existance of)
#' the template (and complement) data.#' 
#' 
#' @importFrom stringr str_detect
#' @keywords internal
#' @noRd 
.fast5status <- function(files, warn = FALSE) {
    
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
    
    ## is /Analysis/BaseCalled_2D present?
    basecalled_2d <- sapply(lsList, FUN = function(x) {
        select(x, name) %>% 
            str_detect(pattern = "BaseCalled_2D")
    })
    
    loc_template <- sapply(lsList, .strandExistence, strand = "BaseCalled_template")
    loc_template <- .chooseStrand(loc_template)
    
    loc_complement <- sapply(lsList, .strandExistence, strand = "BaseCalled_complement")
    loc_complement <- .chooseStrand(loc_complement)
    
    loc_2D <- sapply(lsList, .strandExistence, strand = "BaseCalled_2D")
    loc_2D <- .chooseStrand(loc_2D)

    return(list(read_in_name = all(readInName),
                raw_reads = all(rawReads),
                event_detection = all(eventDetection),
                loc_template = loc_template,
                loc_complement = loc_complement,
                loc_2D = loc_2D,
                basecalled_2d = any(basecalled_2d))) 
    
}

