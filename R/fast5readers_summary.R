.getRawSummary <- function(file, readNumber = NULL, dontCheck = FALSE) {
    
    if(is.character(file)) {
        fid <- H5Fopen(file)
        on.exit(H5Fclose(fid))
    } else {
        fid <- file
    }
    
    if(dontCheck || .groupExistsObj(fid, group = "/Raw/Reads/")) {
        ## get the Read_No. if we need to
        if(is.null(readNumber)) {
            readNumber <- .getReadNumber(fid)
        }
        
        ## Open the group and read the two attributes we want
        gid <- H5Gopen(fid, paste0("/Raw/Reads/Read_", readNumber)) 
        #aid <- H5Aopen(gid, "duration")
        #duration <- H5Aread(aid) 
        #H5Aclose(aid)
        #aid <- H5Aopen(gid, "start_time") 
        #start_time <- H5Aread(aid) 
        #H5Aclose(aid)   
        
        did <- H5Dopen(gid, "Signal")
        signal <- as.integer(H5Dread(did, bit64conversion = "int", compoundAsDataFrame = FALSE))
        median_signal <- median(signal)
        duration <- length(signal)
        H5Dclose(did)
        
        H5Gclose(gid)
    } else {
        median_signal <- duration <- NA
    } 
    
    return(tibble(median_signal, duration))  
}

.getEventsSummary <- function(file, readNumber = NULL, dontCheck = FALSE) {
    
    if(is.character(file)) {
        fid <- H5Fopen(file)
        on.exit(H5Fclose(fid))
    } else {
        fid <- file
    }

    if(dontCheck || .groupExistsObj(fid, group = "/Analyses/EventDetection_000/Reads")) {
        ## get the Read_No. if we need to
        if(is.null(readNumber)) {
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

.getBaseCalledSummary <- function(file, strand = "template") {
    
    if(is.character(file)) {
        fid <- H5Fopen(file)
        on.exit(H5Fclose(fid))
    } else {
        fid <- file
    }
    
    analysisNum <- .findAnalysisNumber(fid)
    
    ## we have to cope with data being under either 1D or 2D file structure
    for(d in c("1D", "2D")) {
        exists <- .groupExistsObj(fid, group = paste0("/Analyses/Basecall_", d, "_", analysisNum, "/Summary/basecall_1d_", strand))
        if(exists) break;
    }
    
    if(!exists) {
        num_events <- duration <- start_time <- NA
    } else {
        ## Open the group and read the attribute we want
        gid <- H5Gopen(fid, paste0("/Analyses/Basecall_", d, "_", analysisNum, "/Summary/basecall_1d_", strand))
        aid <- H5Aopen(gid, "num_events")
        num_events <- H5Aread(aid)
        H5Aclose(aid)   
        H5Gclose(gid)
        
        did <- H5Dopen(fid, paste0("/Analyses/Basecall_", d, "_", analysisNum, "/BaseCalled_", strand, "/Events"))   
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

