## Check whether group is present in a fast5 file
## 
## file Character string containing the name of the file to be read.
## group Character string holding the group we are checking for.
.groupExistsString <- function(file, group) {
    fid <- H5Fopen(file)
    res <- H5Lexists(fid, group)
    H5Fclose(fid)
    return(res)
}

## Check whether group is present in a fast5 file
## 
## @param fid Object of class H5IdComponent created with H5Fopen().
## @param group Character string holding the group we are checking for.
.groupExistsObj <- function(fid, group) {
    res <- H5Lexists(fid, group)
    return(res)
}

## Some fast5 files can be corrupt.  This just checks whether we can open them
## and returns a boolean result.
#' @importFrom methods is
.checkOpening <- function(file) {
    fid <- .H5Fopen_tryCatch(file)
    res <- is(fid,"H5IdComponent")
    if(res) { H5Fclose(fid) }
    return(res)
}

.H5Fopen_tryCatch <- function(file) {
    out <- tryCatch(
        {
            H5Fopen(file)
        },
        error=function(cond) {
            message(paste('Error opening:', file, "- File Skipped"))
            return(NA)
        }
    )    
    return(out)
}

## Determines whether a read is in a pass or fail folder by looking at the path.
## If neither term is found NA is returned
.passFailStatus <- function(path) {
    if(!grepl('pass|fail', path)) {
        return(NA)
    } else {
        return( grepl('pass', path) )
    }
}

.getRunID <- function(file) {
    
    fid <- H5Fopen(file)
    on.exit(H5Fclose(fid))
    
    exists <- .groupExistsObj(fid, group = "/UniqueGlobalKey/tracking_id")
    if(!exists) {
        run_id <- NA
    } else {
        gid <- H5Gopen(fid, "/UniqueGlobalKey/tracking_id")   
        aid <- H5Aopen(gid, "run_id")
        run_id <- H5Aread(aid)
        H5Aclose(aid)
        H5Gclose(gid)
    }
    return( run_id )
    
}

## The sampling is how many times the signal is recorded per second.
## We use this to convert the 'duration' and 'start_time' in the event data
## into seconds.  It may also be useful meta data
#' @importFrom stats median
.getSamplingRate <- function(file) {
    
    fid <- H5Fopen(file)
    on.exit(H5Fclose(fid))
    
    exists <- .groupExistsObj(fid, group = "/UniqueGlobalKey/channel_id")
    if(!exists) {
        rate <- NA
    } else {
        gid <- H5Gopen(fid, "/UniqueGlobalKey/channel_id/")   
        aid <- H5Aopen(gid, "sampling_rate")
        rate <- H5Aread(aid)
        H5Aclose(aid)
        H5Gclose(gid)
    }
    return( rate )
}


.getReadChannelMux <- function(file) {
    
    if(is.character(file)) {
        fid <- H5Fopen(file)
        on.exit(H5Fclose(fid))
    } else {
        fid <- file
    }
    
    # here we get the starting mux and the read number from the channel
    exists <- .groupExistsObj(fid, group = "/Analyses/EventDetection_000/Reads")
    if(!exists) {
        start_mux <- NA
        read_number <- NA
    } else {
        ## get the Read_No., this changes in every file
        gid <- H5Gopen(fid, "/Analyses/EventDetection_000/Reads")
        read_number_char <- h5ls(gid)[1,"name"]
        H5Gclose(gid)
        
        ## Open the group and read the two attribute we want
        gid <- H5Gopen(fid, paste0("/Analyses/EventDetection_000/Reads/", read_number_char))   
        aid <- H5Aopen(gid, "start_mux")
        start_mux <- H5Aread(aid)
        H5Aclose(aid)
        aid <- H5Aopen(gid, "read_number")
        read_number <- H5Aread(aid)
        H5Aclose(aid)   
        H5Gclose(gid)
    }
    
    # we're also interested in the channel information
    exists <- .groupExistsObj(fid, group = "/UniqueGlobalKey/channel_id")
    if(!exists) {
        channel_number <- NA
    } else {
        gid <- H5Gopen(fid, "/UniqueGlobalKey/channel_id/")   
        aid <- H5Aopen(gid, "channel_number")
        channel_number <- H5Aread(aid)
        H5Aclose(aid)
        H5Gclose(gid)
    }
    
    return( tibble(read = as.integer(read_number), channel = as.integer(channel_number), mux = as.integer(start_mux)) ) 
}

.getSummaryEvents <- function(file) {
    
    if(is.character(file)) {
        fid <- H5Fopen(file)
        on.exit(H5Fclose(fid))
    } else {
        fid <- file
    }
    
    exists <- .groupExistsObj(fid, group = "/Analyses/EventDetection_000/Reads")
    if(!exists) {
        start_time <- duration <- num_events <- median_signal <- NA
    } else {
        ## get the Read_No., this changes in every file
        gid <- H5Gopen(fid, "/Analyses/EventDetection_000/Reads")
        read_number_char <- h5ls(gid)[1,"name"]
        H5Gclose(gid)
        
        ## Open the group and read the two attributes we want
        gid <- H5Gopen(fid, paste0("/Analyses/EventDetection_000/Reads/", read_number_char)) 
        aid <- H5Aopen(gid, "duration")
        duration <- H5Aread(aid) 
        H5Aclose(aid)
        aid <- H5Aopen(gid, "start_time") 
        start_time <- H5Aread(aid) 
        H5Aclose(aid)   
        
        did <- H5Dopen(gid, "Events")
        events <- H5Dread(did, bit64conversion = "int", compoundAsDataFrame = FALSE)$mean
        median_signal <- median(events)
        num_events <- length(events)
        H5Dclose(did)
        
        H5Gclose(gid)
    } 
    
    return(data.frame(start_time, duration, num_events, median_signal))  
}

.getEvents <- function(file) {
    
    fid <- H5Fopen(file)
    on.exit(H5Fclose(fid))
    
    exists <- .groupExistsObj(fid, group = "/Analyses/EventDetection_000/Reads")
    if(!exists) {
        start_time <- duration <- num_events <- median_signal <- NA
    } else {
        ## get the Read_No., this changes in every file
        gid <- H5Gopen(fid, "/Analyses/EventDetection_000/Reads")
        read_number_char <- h5ls(gid)[1,"name"]
        H5Gclose(gid)
        
        ## Open the group
        gid <- H5Gopen(fid, paste0("/Analyses/EventDetection_000/Reads/", read_number_char)) 
        did <- H5Dopen(gid, "Events")
        events <- as_tibble(H5Dread(did, bit64conversion = "int", compoundAsDataFrame = TRUE))
        H5Dclose(did)
        
        H5Gclose(gid)
    } 
    
    return(events)
}

.getSummaryBaseCalled <- function(file, strand = "template") {
    
    if(is.character(file)) {
        fid <- H5Fopen(file)
        on.exit(H5Fclose(fid))
    } else {
        fid <- file
    }
    
    ## we have to cope with data being under either 1D or 2D file structure
    for(d in c("1D", "2D")) {
        exists <- .groupExistsObj(fid, group = paste0("/Analyses/Basecall_", d, "_000/Summary/basecall_1d_", strand))
        if(exists) break;
    }
    
    if(!exists) {
        num_events <- duration <- start_time <- NA
    } else {
        ## Open the group and read the attribute we want
        gid <- H5Gopen(fid, paste0("/Analyses/Basecall_", d, "_000/Summary/basecall_1d_", strand))
        aid <- H5Aopen(gid, "num_events")
        num_events <- H5Aread(aid)
        H5Aclose(aid)   
        H5Gclose(gid)
        
        did <- H5Dopen(fid, paste0("/Analyses/Basecall_", d, "_000/BaseCalled_", strand, "/Events"))   
        aid <- H5Aopen(did, "duration")
        duration <- H5Aread(aid)
        H5Aclose(aid)
        aid <- H5Aopen(did, "start_time")
        start_time <- H5Aread(aid)
        H5Aclose(aid)   
        H5Dclose(did)
    }
    basecalledStats <- tibble(num_events, duration, start_time, strand)
    
    return(basecalledStats)  
}

.getBaseCalled <- function(file, strand = "template") {
    
    fid <- H5Fopen(file)
    on.exit(H5Fclose(fid))
    
    for(d in c("1D", "2D")) {
        exists <- .groupExistsObj(fid, group = paste0("/Analyses/Basecall_", d, "_000/Summary/basecall_1d_", strand))
        if(exists) break;
    }
    
    if(!exists) {
        events <- NA
    } else {
        ## Open the group and read
        did <- H5Dopen(fid, paste0("/Analyses/Basecall_", d, "_000/BaseCalled_", strand, "/Events"))   
        events <- as_tibble(H5Dread(did, bit64conversion = "int", compoundAsDataFrame = TRUE))
        H5Dclose(did)
    }

    return(events)  
}

.getAligned <- function(file) {
 
    fid <- H5Fopen(file)
    on.exit(H5Fclose(fid))
    
    exists <- .groupExistsObj(fid, group = "/Analyses/Alignment_000/Aligned_2d/SAM")
    
    if(exists) {
        did <- H5Dopen(fid, "/Analyses/Alignment_000/Aligned_2d/SAM")
        samData <- data.table(H5Dread(did))
    } else {
        samData <- ""
    }
    
    return(samData)
}



#' Read the log entry from a single fast5 file
#' 
#' Basecalling procedures performed on fast5 files generally leave a 
#' log file entry recording how far through the pipelin the file 
#' proceeded.  This function will extract this information as a 
#' single string.  It can be printed in a more readable format 
#' using the \link{\code{cat()}} function.
#' 
#' @param file Character vector of fast5 file to be read.
#' @return Character vector containing the log information.  
#' \code{NULL} if no log is found.
#' @examples 
#' fast5file <- system.file('extdata', c('example.fast5', package = "IONiseR")
#' log <- readFast5Log(fast5file)
#' cat(tmp)
#' @export
readFast5Log <- function(file) {

    fid <- H5Fopen(file)
    on.exit(H5Fclose(fid))
    
    log <- NULL
    
    for(d in c("1D", "2D")) {
        path <- paste0("/Analyses/Basecall_", d, "_000/Log")
        exists <- .groupExistsObj(fid, group = path)
        if(exists) {
            did <- H5Dopen(fid, path)
            log <- c(log, H5Dread(did))
        }
    }
    
    if(is.null(log))
        message("No log information found.")
    
    return(log)
}

