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
    
    ## some files can't be opened, so we filter them here
    message("Checking file validity")
    fileStatus <- sapply(files, .checkOpening, USE.NAMES = FALSE)
    files <- files[ which(fileStatus) ]
    
    message("Reading Channel Data")
    readInfo <- lapply(files, .getReadChannelMux)
    readInfo <- rbindlist(readInfo)
    readInfo <- data.table(id = 1:nrow(readInfo), file = basename(files), readInfo)
    
    message("Reading Raw Data")
    rawData <- lapply(files, .getSummaryRaw)
    rawData <- data.table(id = readInfo[,id], rbindlist(rawData))
    ## we convert timing data into seconds. 
    ## To do this we find the sampling rate stored in one file
    samplingRate <- .getSamplingRate(files[1])
    rawData <- mutate(rawData, 
                      start_time = start_time / samplingRate,
                      duration = duration / samplingRate)
    
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
    fq_t <- sapply(files[ template[,'id'] ], .getFastqString, strand = "template")
    fq_t <- .processFastqVec(fq_t, readIDs = template[, 'id'], appendID = "_template")  
    fastq_template <- fq_t$fastq
    ## if there are any invalid entries we need to remove them
    if(length(fq_t$invalid)) {
        template <- template[-fq_t$invalid,]
    }
    
    if(nrow(complement)) {
        message("Reading Complement FASTQ")
        fq_c <- sapply(files[ complement[,'id'] ], .getFastqString, strand = "complement")
        fq_c <- .processFastqVec(fq_c, readIDs = complement[, 'id'], appendID = "_complement")  
        fastq_complement <- fq_c$fastq
        ## if there are any invalid entries we need to remove them
        if(length(fq_c$invalid)) {
            complement <- complement[-fq_c$invalid,]
        }
    } else {
        fastq_complement <- ShortRead::ShortReadQ()
    }
    
    ## we haven't read anything about 2D reads yet, so we need to identify which
    ## files have them.  Then we'll read only those
    idx2D <- which(sapply(files, .groupExistsString, group = paste0("/Analyses/Basecall_2D_000/BaseCalled_2D")))
    if(length(idx2D)) {
        message("Reading 2D FASTQ")
        fq_2D <- sapply(files[ idx2D ], .getFastqString, strand = "2D")
        fq_2D <- .processFastqVec(fq_2D, readIDs = readInfo[idx2D, id], appendID = "_2D")  
        fastq_2D <- fq_2D$fastq
    }

    ## We update the individual strands to indicate if they are part of a full 2D read
    template <- data.table(template, full_2D = template[,'id'] %in% idx2D)
    if(nrow(complement)) {
        complement <- data.table(complement, full_2D = complement[,'id'] %in% idx2D)
    }
    
    ## combine the template, complement and 2D data
    baseCalled <- rbind(template, complement)
    fastq <- append(fastq_template, fastq_complement)
    if(length(idx2D)) {
        fastq <- append(fastq, fastq_2D)
    }
    
    message("Done")
    obj <- new("Fast5Summary", readInfo = readInfo, rawData = data.table(rawData), baseCalled = baseCalled, fastq = fastq)
    
    return(obj)
}



## work file by file, rather than by data type
## doens't seem to be any faster though
readFast5Summary2 <- function(files) {
    
    readInfo <- data.table(id = integer(length = length(files)),
                           file = character(length = length(files)),
                           read = integer(length = length(files)), 
                           channel = integer(length = length(files)),
                           mux = integer(length = length(files)),
                           pass = logical(length = length(files)))
    
    rawData <- data.table(id = integer(length = length(files)),
                          start_time = numeric(length = length(files)),
                          duration = numeric(length = length(files)), 
                          num_events = integer(length = length(files)),
                          median_signal = numeric(length = length(files)))
    
    templateData <- complementData <- 
        data.table(id = integer(length = length(files)),
                   num_events = integer(length = length(files)),
                   duration = numeric(length = length(files)), 
                   start_time = numeric(length = length(files)),
                   strand = character(length = length(files)))
    
    fastq_t_vec <- fastq_c_vec <- fastq_2d_vec <- character(length = length(files))
            
    samplingRate <- .getSamplingRate(files[1])
    
    pb = txtProgressBar(min = 0, max = length(files), initial = 0, style = 3) 
    
    for(i in 1:length(files)) {
        
        setTxtProgressBar(pb, value = i)
        
        fileStatus <- suppressMessages(.checkOpening( files[i] ))
        if(!fileStatus) {
            next()
        }
        
        fid <- H5Fopen(files[i])

        ## Read Channel Data
        rcm <- .getReadChannelMux( fid )
        readInfo[i,] <- cbind(i, basename(files[i]), rcm, .passFailStatus(files[i]))
        ## Read Raw Data
        raw <- .getSummaryRaw( fid )
        rawData[i,] <- cbind(i, raw[1:2] / samplingRate, raw[3:4])
        
        ## Template data
        template <- .getSummaryBaseCalled(fid, strand = "template")
        templateData[i,] <- cbind(i, template)
        ## Complement data
        complement <- .getSummaryBaseCalled(fid, strand = "complement")
        complementData[i,] <- cbind(i, complement)
        
        if(!is.na(template[1,num_events]))
            fastq_t_vec[i] <- .getFastqString(fid, strand = "template")
        if(!is.na(complement[1,num_events]))
            fastq_c_vec[i] <- .getFastqString(fid, strand = "complement")
        
        ## reading 2D fastq
        if(.groupExistsObj(fid, group = paste0("/Analyses/Basecall_2D_000/BaseCalled_2D"))) {
            fastq_2d_vec[i] <- .getFastqString(fid, strand = "2D")
        }
        H5Fclose(fid)
    }
    ## if we skipped any files, tell the user
    if(length(which(readInfo[,id] == 0))) {
        message("Error open the following files, files were ignored:")
        message(paste0(files[ which(readInfo[,id] == 0) ], collapse = "\n"))
    }
    ## filter entries where no data were recorded
    readInfo <- filter(readInfo, id > 0)
    rawData <- filter(rawData, id > 0)
    templateData <- filter(templateData, !is.na(start_time), id > 0)
    complementData <- filter(complementData, !is.na(start_time), id > 0)
    ## We update the individual strands to indicate if they are part of a full 2D read
    templateData <- mutate(templateData, full_2D = templateData[,id] %in% which(nchar(fastq_2d_vec) > 0))
    complementData <- mutate(complementData, full_2D = complementData[,id] %in% which(nchar(fastq_2d_vec) > 0))
    
    ## process the fastq strings into objects
    fq_t <- .processFastqVec(fastq_t_vec, readIDs = templateData[,id], appendID = "_template")
    fastq_template <- fq_t$fastq
    fq_c <- .processFastqVec(fastq_c_vec, readIDs = complementData[,id], appendID = "_complement")
    fastq_complement <- fq_c$fastq
    fq_2D <- .processFastqVec(fastq_2d_vec, readIDs = which(nchar(fastq_2d_vec) > 0), appendID = "_2D")
    fastq_2D <- fq_2D$fastq
    ## combine the template, complement and 2D data
    baseCalled <- rbind(templateData, complementData)
    fastq <- append(fastq_template, fastq_complement)
    fastq <- append(fastq, fastq_2D)
    
    obj <- new("Fast5Summary", readInfo = readInfo, rawData = rawData, baseCalled = baseCalled, fastq = fastq)
    
    return(obj)
}


## work file by file, rather than by data type
#' @importFrom utils setTxtProgressBar txtProgressBar
readFast5Summary3 <- function(files) {
    
    readInfo <- data.table(id = integer(length = length(files)),
                           file = character(length = length(files)),
                           read = integer(length = length(files)), 
                           channel = integer(length = length(files)),
                           mux = integer(length = length(files)),
                           pass = logical(length = length(files)))
    
    rawData <- data.table(id = integer(length = length(files)),
                          start_time = numeric(length = length(files)),
                          duration = numeric(length = length(files)), 
                          num_events = integer(length = length(files)),
                          median_signal = numeric(length = length(files)))
    
    templateData <- complementData <- 
        data.table(id = integer(length = length(files)),
                   num_events = integer(length = length(files)),
                   duration = numeric(length = length(files)), 
                   start_time = numeric(length = length(files)),
                   strand = character(length = length(files)))
    
    fastq_t_vec <- fastq_c_vec <- fastq_2d_vec <- character(length = length(files))
    
    samplingRate <- .getSamplingRate(files[1])
    
    pb = txtProgressBar(min = 0, max = length(files), initial = 0, style = 3) 
    
    for(i in 1:length(files)) {
        
        setTxtProgressBar(pb, value = i)
        
        fileStatus <- suppressMessages(.checkOpening( files[i] ))
        if(!fileStatus) {
            next()
        }
        
        fid <- H5Fopen(files[i])
        
        ## Read Channel Data
        rcm <- .getReadChannelMux( fid )
        readInfo[i,] <- cbind(i, basename(files[i]), rcm, .passFailStatus(files[i]))
        ## Read Raw Data
        raw <- .getSummaryRaw( fid )
        rawData[i,1] <- i
        rawData[i,2:3] <- raw[i,1:2] / samplingRate
        rawData[i,4:5] <- raw[3:4] 
        
        ## Template data
        template <- .getSummaryBaseCalled(fid, strand = "template")
        templateData[i,] <- cbind(i, template)
        ## Complement data
        complement <- .getSummaryBaseCalled(fid, strand = "complement")
        complementData[i,1] <- i
        complementData[i,2:5] <- complement
        
        if(!is.na(template[1,num_events]))
            fastq_t_vec[i] <- .getFastqString(fid, strand = "template")
        if(!is.na(complement[1,num_events]))
            fastq_c_vec[i] <- .getFastqString(fid, strand = "complement")
        
        ## reading 2D fastq
        if(.groupExistsObj(fid, group = paste0("/Analyses/Basecall_2D_000/BaseCalled_2D"))) {
            fastq_2d_vec[i] <- .getFastqString(fid, strand = "2D")
        }
        H5Fclose(fid)
    }
    ## if we skipped any files, tell the user
    if(length(which(readInfo[,id] == 0))) {
        message("Error open the following files, files were ignored:")
        message(paste0(files[ which(readInfo[,id] == 0) ]), collapse = "\n")
    }
    ## filter entries where no data were recorded
    readInfo <- filter(readInfo, id > 0)
    rawData <- filter(rawData, id > 0)
    templateData <- filter(templateData, !is.na(start_time), id > 0)
    complementData <- filter(complementData, !is.na(start_time), id > 0)
    ## We update the individual strands to indicate if they are part of a full 2D read
    templateData <- mutate(templateData, full_2D = templateData[,id] %in% which(nchar(fastq_2d_vec) > 0))
    complementData <- mutate(complementData, full_2D = complementData[,id] %in% which(nchar(fastq_2d_vec) > 0))
    
    ## process the fastq strings into objects
    fq_t <- .processFastqVec(fastq_t_vec, readIDs = templateData[,id], appendID = "_template")
    fastq_template <- fq_t$fastq
    fq_c <- .processFastqVec(fastq_c_vec, readIDs = complementData[,id], appendID = "_complement")
    fastq_complement <- fq_c$fastq
    fq_2D <- .processFastqVec(fastq_2d_vec, readIDs = which(nchar(fastq_2d_vec) > 0), appendID = "_2D")
    fastq_2D <- fq_2D$fastq
    ## combine the template, complement and 2D data
    baseCalled <- rbind(templateData, complementData)
    fastq <- append(fastq_template, fastq_complement)
    fastq <- append(fastq, fastq_2D)
    
    obj <- new("Fast5Summary", readInfo = readInfo, rawData = rawData, baseCalled = baseCalled, fastq = fastq)
    
    return(obj)
}



## work file by file, rather than by data type
#' @importFrom utils untar
readFast5SummaryTar <- function(tarfile, chunkSize = 100) {
    
    ## find the path to any fast5 files in the tar archive
    message("Identifying fast5 files in archive")
    files <- untar(tarfile = tarfile, list = TRUE)
    files <- grep(".fast5$", files, value = TRUE)
    ## do we have any files starting with a '.'? Remove them if so
    hidden <- grep("^\\.",basename(files))
    if(length(hidden)) {
        files <- files[ -hidden ]
    }
    
    ## predefine the structures to hold the data
    readInfo <- data.table(id = integer(length = length(files)),
                           file = character(length = length(files)),
                           read = integer(length = length(files)), 
                           channel = integer(length = length(files)),
                           mux = integer(length = length(files)),
                           pass = logical(length = length(files)))
    
    rawData <- data.table(id = integer(length = length(files)),
                          start_time = numeric(length = length(files)),
                          duration = numeric(length = length(files)), 
                          num_events = integer(length = length(files)),
                          median_signal = numeric(length = length(files)))
    
    templateData <- complementData <- 
        data.table(id = integer(length = length(files)),
                   num_events = integer(length = length(files)),
                   duration = numeric(length = length(files)), 
                   start_time = numeric(length = length(files)),
                   strand = character(length = length(files)))
    
    fastq_t_vec <- fastq_c_vec <- fastq_2d_vec <- character(length = length(files))
    
    td <- tempdir()
    pb = txtProgressBar(min = 0, max = length(files), initial = 0, style = 3) 
    
    ## chunk counter
    j <- 1
    while(j <= length(files)) {
        ## untar a block of files, making sure we dont try too many
        end <- min(j+chunkSize-1, length(files))
        untar(tarfile = tarfile, files = files[j:end], exdir = td)
        chunk_files <- file.path(td, files[j:end])
        samplingRate <- .getSamplingRate(chunk_files[1])
        for(i in 1:length(chunk_files)) {
            
            idx <- j + i - 1
            
            fileStatus <- suppressMessages(.checkOpening( chunk_files[i] ))
            if(!fileStatus) {
                next()
            }
            
            fid <- H5Fopen(chunk_files[i])
            
            ## Read Channel Data
            rcm <- .getReadChannelMux( fid )
            readInfo[idx,] <- cbind(idx, basename(chunk_files[i]), rcm, .passFailStatus(chunk_files[i]))
            ## Read Raw Data
            raw <- .getSummaryRaw( fid )
            rawData[idx,] <- cbind(idx, raw[1:2] / samplingRate, raw[3:4])
            
            ## Template data
            template <- .getSummaryBaseCalled(fid, strand = "template")
            templateData[idx,] <- cbind(idx, template)
            ## Complement data
            complement <- .getSummaryBaseCalled(fid, strand = "complement")
            complementData[idx,] <- cbind(idx, complement)
            
            if(!is.na(template[1,num_events]))
                fastq_t_vec[idx] <- .getFastqString(fid, strand = "template")
            if(!is.na(complement[1,num_events]))
                fastq_c_vec[idx] <- .getFastqString(fid, strand = "complement")
            
            ## reading 2D fastq
            if(.groupExistsObj(fid, group = paste0("/Analyses/Basecall_2D_000/BaseCalled_2D"))) {
                fastq_2d_vec[idx] <- .getFastqString(fid, strand = "2D")
            }
            H5Fclose(fid)
        }
        
        
        ## increment our chunk counter and remove temp files
        removeStatus <- file.remove(file.path(td, files[j:end]))
        j <- end + 1
        setTxtProgressBar(pb, value = j-1)
    }

    ## if we skipped any files, tell the user
    if(length(which(readInfo[,id] == 0))) {
        message("Error open the following files, files were ignored:")
        message(paste0(files[ which(readInfo[,id] == 0) ]), collapse = "\n")
    }
    ## filter entries where no data were recorded
    readInfo <- filter(readInfo, id > 0)
    rawData <- filter(rawData, id > 0)
    templateData <- filter(templateData, !is.na(start_time), id > 0)
    complementData <- filter(complementData, !is.na(start_time), id > 0)
    ## We update the individual strands to indicate if they are part of a full 2D read
    templateData <- mutate(templateData, full_2D = templateData[,id] %in% which(nchar(fastq_2d_vec) > 0))
    complementData <- mutate(complementData, full_2D = complementData[,id] %in% which(nchar(fastq_2d_vec) > 0))
    
    ## process the fastq strings into objects
    fq_t <- .processFastqVec(fastq_t_vec, readIDs = templateData[,id], appendID = "_template")
    fastq_template <- fq_t$fastq
    fq_c <- .processFastqVec(fastq_c_vec, readIDs = complementData[,id], appendID = "_complement")
    fastq_complement <- fq_c$fastq
    fq_2D <- .processFastqVec(fastq_2d_vec, readIDs = which(nchar(fastq_2d_vec) > 0), appendID = "_2D")
    fastq_2D <- fq_2D$fastq
    ## combine the template, complement and 2D data
    baseCalled <- rbind(templateData, complementData)
    fastq <- append(fastq_template, fastq_complement)
    fastq <- append(fastq, fastq_2D)
    
    obj <- new("Fast5Summary", readInfo = readInfo, rawData = rawData, baseCalled = baseCalled, fastq = fastq)
    
    return(obj)

}