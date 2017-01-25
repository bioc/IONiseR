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
        signal <- as.integer(H5Dread(did, bit64conversion = "int", compoundAsDataFrame = FALSE))
        median_signal <- median(signal)
        duration <- length(signal)
        H5Dclose(did)
        
        #H5Gclose(gid)
    } else {
        median_signal <- duration <- NA
    } 
    
    return(tibble(median_signal, duration))  
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
        aid <- H5Aopen(gid, "start_time") 
        start_time <- H5Aread(aid) 
        H5Aclose(aid)
        H5Gclose(gid)
        
        gid <- H5Gopen(fid, paste0("/Analyses/EventDetection_000/Summary/event_detection")) 
        aid <- H5Aopen(gid, "num_events") 
        num_events <- H5Aread(aid) 
        H5Aclose(aid)
        H5Gclose(gid)
    } else {
        start_time <- num_events <- NA
    }
    
    return(tibble(start_time, num_events))  
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
        
        aid <- H5Aopen(gid, "num_events")
        num_events <- H5Aread(aid)
        H5Aclose(aid)   
        aid <- H5Aopen(gid, "num_skips")
        num_skips <- H5Aread(aid)
        H5Aclose(aid)  
        aid <- H5Aopen(gid, "num_stays")
        num_stays <- H5Aread(aid)
        H5Aclose(aid)  
        
        H5Gclose(gid)
    }
    basecalledStats <- tibble(num_events, num_skips, num_stays, strand)
    
    return(basecalledStats)  
}

