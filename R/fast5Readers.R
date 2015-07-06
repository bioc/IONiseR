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

## The sampling is how many times the signal is recorded per second.
## We use this to convert the 'duration' and 'start_time' in the raw data
## into seconds.  It may also be useful meta data
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
    
    fid <- H5Fopen(file)
    on.exit(H5Fclose(fid))
    
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
    
    return( data.frame(read = as.integer(read_number), channel = as.integer(channel_number), mux = as.integer(start_mux)) ) 
}

.getSummaryRaw <- function(file) {
    
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

.getRaw <- function(file) {
    
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
        events <- data.table(H5Dread(did, bit64conversion = "int", compoundAsDataFrame = TRUE))
        H5Dclose(did)
        
        H5Gclose(gid)
    } 
    
    return(events)
}

.getSummaryBaseCalled <- function(file, strand = "template") {
    
    fid <- H5Fopen(file)
    on.exit(H5Fclose(fid))
    
    exists <- .groupExistsObj(fid, group = paste0("/Analyses/Basecall_2D_000/Summary/basecall_1d_", strand))
    if(!exists) {
        num_events <- duration <- start_time <- NA
    } else {
        ## Open the group and read the attribute we want
        gid <- H5Gopen(fid, paste0("/Analyses/Basecall_2D_000/Summary/basecall_1d_", strand))
        aid <- H5Aopen(gid, "num_events")
        num_events <- H5Aread(aid)
        H5Aclose(aid)   
        H5Gclose(gid)
        
        did <- H5Dopen(fid, paste0("/Analyses/Basecall_2D_000/BaseCalled_", strand, "/Events"))   
        aid <- H5Aopen(did, "duration")
        duration <- H5Aread(aid)
        H5Aclose(aid)
        aid <- H5Aopen(did, "start_time")
        start_time <- H5Aread(aid)
        H5Aclose(aid)   
        H5Dclose(did)
    }
    basecalledStats <- data.table(num_events, duration, start_time, strand)
    
    return(basecalledStats)  
}

.getBaseCalled <- function(file, strand = "template") {
    
    fid <- H5Fopen(file)
    on.exit(H5Fclose(fid))
    
    exists <- .groupExistsObj(fid, group = paste0("/Analyses/Basecall_2D_000/Summary/basecall_1d_", strand))
    if(!exists) {
        events <- NA
    } else {
        ## Open the group and read
        did <- H5Dopen(fid, paste0("/Analyses/Basecall_2D_000/BaseCalled_", strand, "/Events"))   
        events <- data.table(H5Dread(did, bit64conversion = "int", compoundAsDataFrame = TRUE))
        H5Dclose(did)
    }

    return(events)  
}