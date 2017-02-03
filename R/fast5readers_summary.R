#' @import bit64
.getRawSummary <- function(file, readNumber = NA, dontCheck = FALSE) {
    
    if(is.character(file)) {
        fid <- H5Fopen(file)
        on.exit(H5Fclose(fid))
    } else {
        fid <- file
    }
    
    if(dontCheck || .groupExistsObj(fid, group = "/Raw/Reads/")) {
        ## get the Read_No. if we need to
        if(is.na(readNumber)) {
            readNumber <- .getReadNumber(fid)
        }
        
        ## Open the group and read the two attributes we want
        #gid <- H5Gopen(fid, paste0("/Raw/Reads/Read_", readNumber)) 
        #aid <- H5Aopen(gid, "duration")
        #duration <- H5Aread(aid) 
        #H5Aclose(aid)
        #aid <- H5Aopen(gid, "start_time") 
        #start_time <- H5Aread(aid) 
        #H5Aclose(aid)   
        
        did <- H5Dopen(fid, paste0("/Raw/Reads/Read_", readNumber, "/Signal"))
        #signal <- H5Dread(did, bit64conversion = "int", compoundAsDataFrame = FALSE)
        signal <- H5Dread(did, compoundAsDataFrame = FALSE)
        H5Dclose(did)
        mean_signal <- mean(signal)
        #duration <- length(signal)

        #H5Gclose(gid)
    } else {
        mean_signal <- NA
    } 
    
    return(tibble(mean_signal))  
}

.getEventsSummary <- function(file, readNumber = NA, dontCheck = FALSE) {
    
    if(is.character(file)) {
        fid <- H5Fopen(file)
        on.exit(H5Fclose(fid))
    } else {
        fid <- file
    }

    if(dontCheck || .groupExistsObj(fid, group = "/Analyses/EventDetection_000/Reads")) {
        ## get the Read_No. if we need to
        if(is.na(readNumber)) {
            readNumber <- .getReadNumber(fid)
        }

        ## Open the group and read the attributes we want
        gid <- H5Gopen(fid, paste0("/Analyses/EventDetection_000/Reads/Read_", readNumber)) 
        start_time <- tryCatch(.readAttribute(gid, attribute = "start_time"), 
                             error = function(e) { NA })
        duration <- tryCatch(.readAttribute(gid, attribute = "duration"), 
                             error = function(e) { NA })
        H5Gclose(gid)
        
        if(.groupExistsObj(fid, paste0("/Analyses/EventDetection_000/Summary/event_detection"))) {
            gid <- H5Gopen(fid, paste0("/Analyses/EventDetection_000/Summary/event_detection")) 
            num_events <- tryCatch(.readAttribute(gid, attribute = "num_events"), 
                                 error = function(e) { NA }) 
            H5Gclose(gid)
        } else {
            did <- H5Dopen(fid, paste0("/Analyses/EventDetection_000/Reads/Read_", readNumber, "/Events"))
            h5dataset <- H5Dget_space(did)    
            num_events <- H5Sget_simple_extent_dims(h5dataset)$size
            H5Dclose(did)
        }
    } else {
        start_time <- duration <- num_events <- NA
    }
    
    return(tibble(start_time, duration, num_events))  
}

.readAttribute <- function(gid, attribute) {
    aid <- H5Aopen(gid, attribute)
    on.exit( H5Aclose(aid) )
    val <- H5Aread(aid)
    return(val)
}


.getBaseCalledSummary <- function(file, strand = "template", d = "1D", analysisNum = "000",
                                  dontCheck = FALSE) {
    
    if(is.character(file)) {
        fid <- H5Fopen(file)
        on.exit(H5Fclose(fid))
    } else {
        fid <- file
    }
    
    if(dontCheck || .groupExistsObj(fid, group = paste0("/Analyses/Basecall_", d, "_", analysisNum, "/Summary/basecall_1d_", strand))) {
        ## Open the group and read the attribute we want
        gid <- H5Gopen(fid, paste0("/Analyses/Basecall_", d, "_", analysisNum, "/Summary/basecall_1d_", strand))
        
        num_events <- tryCatch(.readAttribute(gid, attribute = "num_events"), 
                               error = function(e) { NA })
        num_skips <- tryCatch(.readAttribute(gid, attribute = "num_skips"), 
                              error = function(e) { NA })
        num_stays <- tryCatch(.readAttribute(gid, attribute = "num_stays"), 
                              error = function(e) { NA })
        
        H5Gclose(gid)
    } else {
        num_events <- num_skips <- num_stays <- NA
    }
    basecalledStats <- tibble(num_events, num_skips, num_stays, strand)
    
    return(basecalledStats)  
}

