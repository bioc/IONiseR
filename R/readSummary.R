#' Read summary data from fast5 files.
#' 
#' Reads one or more fast5 files and collects summary information about them.  
#' 
#' Currently this function assumes all files passed to it come from the same 
#' sequencing run.  It makes no effort to check for alternative file names or 
#' the like.  If files from multiple runs are passed to it they will be 
#' collated together and any analysis performed on them will represent the 
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
#' @importFrom utils packageVersion
readFast5Summary <- function(files) {
    
    ## some files can't be opened, so we filter them here
    message("Checking file validity")
    files <- files[file.exists(files)]
    if(length(files) == 0) {
        stop('None of the provided files can be accessed.  ',
             'Have you supplied the correct path?')
    }
    fileStatus <- sapply(files, .checkOpening, USE.NAMES = FALSE)
    files <- files[ which(fileStatus) ]
    
    status <- .fast5status(files = sample(files, size = min(length(files), 15)))
    if(status$template_loc == "") {
        stop("No basecalls for template strand found.  Aborting")
    }
    
    message("Reading Channel Data")
    if(status$read_in_name) {
        readNums <- as.integer(str_match(string = files, pattern = "read([0-9]+)")[,2])
    } else { 
        readNums <- rep(NA, length(files))
    }
    readInfo <- do.call("rbind", mapply(.getReadChannelMux2, files, readNums, dontCheck = TRUE,
                                        USE.NAMES = FALSE, SIMPLIFY = FALSE))
    readInfo <- as_tibble(cbind(id = 1:nrow(readInfo), file = basename(files), readInfo))
    
    if(status$raw_reads) {
        message("Reading Raw Data")
        ## raw data is the median signal, and the duration of the read
        rawData <- as_tibble(do.call("rbind", mapply(.getRawSummary, files, readNums, dontCheck = TRUE, 
                                                     USE.NAMES = FALSE, SIMPLIFY = FALSE)))
    } 
    
    message("Reading Event Data")
    eventData <- do.call("rbind", mapply(.getEventsSummary,  files, readNums, dontCheck = TRUE, 
                                         USE.NAMES = FALSE, SIMPLIFY = FALSE))
    
    if(status$raw_reads) {
        rawEventData <- as_tibble(cbind(id = readInfo[['id']], rawData, eventData))
    } else {
        rawEventData <- as_tibble(cbind(id = readInfo[['id']], eventData))
    }
    
    ## we convert timing data into seconds. 
    ## To do this we find the sampling rate stored in one file
    samplingRate <- .getSamplingRate(files[1])
    rawEventData <- mutate(rawEventData, 
                           start_time = start_time / samplingRate,
                           duration = duration / samplingRate)
    
    message("Reading Template Data")
    d <- str_match(pattern = "_([12]D)_", string = status$template_loc)[,2]
    template <- do.call("rbind", mapply(.getBaseCalledSummary, files, dontCheck = FALSE, 
                                        strand = "template", d = d,
                                        SIMPLIFY = FALSE, USE.NAMES = FALSE))
    template <- mutate(template, id = readInfo[['id']])
    template <- filter(template, !(is.na(num_events)))
    
    message("Reading Template FASTQ")
    ## get the fastq for those that have it
    fq_t <- mapply(.getFastqString, files[ template[['id']] ], strand = "template", d = d)
    fq_t <- .processFastqVec(fq_t, readIDs = template[['id']], appendID = "_template")  
    fastq_template <- fq_t$fastq
    ## if there are any invalid entries we need to remove them
    if(length(fq_t$invalid)) {
        template <- template[-fq_t$invalid,]
    }
    baseCalled <- template
    
    if(nchar(status$complement_loc)) { 
        message("Reading Complement Data")
        d <- str_match(pattern = "_([12]D)_", string = status$complement_loc)[,2]
        complement <- do.call("rbind", mapply(.getBaseCalledSummary, files, dontCheck = FALSE, 
                                              strand = "complement", d = d,
                                              SIMPLIFY = FALSE, USE.NAMES = FALSE))
        complement <- mutate(complement, id = readInfo[['id']])
        complement <- filter(complement, !(is.na(num_events)))
        
        
        message("Reading Complement FASTQ")
        fq_c <- mapply(.getFastqString, files[ complement[['id']] ], strand = "complement", d = d)
        fq_c <- .processFastqVec(fq_c, readIDs = complement[['id']], appendID = "_complement")  
        fastq_complement <- fq_c$fastq
        
        ## if there are any invalid entries we need to remove them
        if(length(fq_c$invalid)) {
            complement <- complement[-fq_c$invalid,]
        }
        
    } else {
        fastq_complement <- ShortRead::ShortReadQ()
        complement <- NULL
    }
    
    ## we haven't read anything about 2D reads yet, so we need to identify which
    ## files have them.  Then we'll read only those
    idx2D <- which(sapply(files, .groupExistsString, group = paste0("/Analyses/Basecall_2D_000/BaseCalled_2D")))
    if(length(idx2D)) {
        message("Reading 2D FASTQ")
        fq_2D <- sapply(files[ idx2D ], .getFastqString, strand = "2D", d = "2D")
        fq_2D <- .processFastqVec(fq_2D, readIDs = readInfo[['id']][idx2D], appendID = "_2D")  
        fastq_2D <- fq_2D$fastq
    }
    
    ## We update the individual strands to indicate if they are part of a full 2D read
    template <- as_tibble(cbind(template, full_2D = template[['id']] %in% idx2D))
    complement <- as_tibble(cbind(complement, full_2D = complement[['id']] %in% idx2D))
    
    # ## combine the template, complement and 2D data
    baseCalled <- rbind(template, complement)
    fastq <- ShortRead::append(fastq_template, fastq_complement)
    if(length(idx2D)) {
        fastq <- ShortRead::append(fastq, fastq_2D)
    }
    
    message("Done")
    obj <- new("Fast5Summary", 
               readInfo = readInfo, 
               rawData = tibble(), 
               eventData = rawEventData, 
               baseCalled = baseCalled, 
               fastq = fastq,
               versions = list('IONiseR' = strsplit(as.character(packageVersion("IONiseR")),".",fixed=TRUE)[[1]]
                               # 'MinKNOW' = max(versions))
               )
    )
    
    return(obj)
}

#' @importFrom BiocParallel MulticoreParam register bpmapply
readFast5Summary.mc <- function(files, ncores = 2) {
  
  register(MulticoreParam(workers = ncores))
  
  ## some files can't be opened, so we filter them here
  message("Checking file validity")
  fileStatus <- sapply(files, .checkOpening, USE.NAMES = FALSE)
  files <- files[ which(fileStatus) ]
  
  status <- .fast5status(files = sample(files, size = min(length(files), 15)))
  if(status$template_loc == "") {
    stop("No basecalling for template strand found.  Aborting")
  }
  
  message("Reading Channel Data")
  if(status$read_in_name) {
    readNums <- as.integer(str_match(string = files, pattern = "read([0-9]+)")[,2])
  } else { 
    readNums <- rep(NA, length(files))
  }
  readInfo <- do.call("rbind", bpmapply(.getReadChannelMux2, files, readNums, dontCheck = TRUE,
                                      USE.NAMES = FALSE, SIMPLIFY = FALSE))
  readInfo <- as_tibble(cbind(id = 1:nrow(readInfo), file = basename(files), readInfo))
  
  if(status$raw_reads) {
    message("Reading Raw Data")
    ## raw data is the median signal, and the duration of the read
    rawData <- as_tibble(do.call("rbind", bpmapply(.getRawSummary, files, readNums, dontCheck = TRUE, 
                                                 USE.NAMES = FALSE, SIMPLIFY = FALSE)))
  } 
  
  message("Reading Event Data")
  eventData <- do.call("rbind", bpmapply(.getEventsSummary,  files, readNums, dontCheck = TRUE, 
                                       USE.NAMES = FALSE, SIMPLIFY = FALSE))
  
  if(status$raw_reads) {
    rawEventData <- as_tibble(cbind(id = readInfo[['id']], rawData, eventData))
  } else {
    rawEventData <- as_tibble(cbind(id = readInfo[['id']], eventData))
  }
  
  ## we convert timing data into seconds. 
  ## To do this we find the sampling rate stored in one file
  samplingRate <- .getSamplingRate(files[1])
  rawEventData <- mutate(rawEventData, 
                         start_time = start_time / samplingRate,
                         duration = duration / samplingRate)
  
  message("Reading Template Data")
  d <- str_match(pattern = "_([12]D)_", string = status$template_loc)[,2]
  template <- do.call("rbind", bpmapply(.getBaseCalledSummary, files, dontCheck = FALSE, 
                                      strand = "template", d = d,
                                      SIMPLIFY = FALSE, USE.NAMES = FALSE))
  template <- mutate(template, id = readInfo[['id']])
  template <- filter(template, !(is.na(num_events)))
  
  message("Reading Template FASTQ")
  ## get the fastq for those that have it
  fq_t <- bpmapply(.getFastqString, files[ template[['id']] ], strand = "template", d = d)
  fq_t <- .processFastqVec(fq_t, readIDs = template[['id']], appendID = "_template")  
  fastq_template <- fq_t$fastq
  ## if there are any invalid entries we need to remove them
  if(length(fq_t$invalid)) {
    template <- template[-fq_t$invalid,]
  }
  baseCalled <- template
  
  if(nchar(status$complement_loc)) { 
    message("Reading Complement Data")
    d <- str_match(pattern = "_([12]D)_", string = status$complement_loc)[,2]
    complement <- do.call("rbind", bpmapply(.getBaseCalledSummary, files, dontCheck = FALSE, 
                                          strand = "complement", d = d,
                                          SIMPLIFY = FALSE, USE.NAMES = FALSE))
    complement <- mutate(complement, id = readInfo[['id']])
    complement <- filter(complement, !(is.na(num_events)))
    
    
    message("Reading Complement FASTQ")
    fq_c <- bpmapply(.getFastqString, files[ complement[['id']] ], strand = "complement", d = d)
    fq_c <- .processFastqVec(fq_c, readIDs = complement[['id']], appendID = "_complement")  
    fastq_complement <- fq_c$fastq
    
    ## if there are any invalid entries we need to remove them
    if(length(fq_c$invalid)) {
      complement <- complement[-fq_c$invalid,]
    }
    
  } else {
    fastq_complement <- ShortRead::ShortReadQ()
    complement <- NULL
  }
  
  ## we haven't read anything about 2D reads yet, so we need to identify which
  ## files have them.  Then we'll read only those
  idx2D <- which(sapply(files, .groupExistsString, group = paste0("/Analyses/Basecall_2D_000/BaseCalled_2D")))
  if(length(idx2D)) {
    message("Reading 2D FASTQ")
    fq_2D <- sapply(files[ idx2D ], .getFastqString, strand = "2D", d = "2D")
    fq_2D <- .processFastqVec(fq_2D, readIDs = readInfo[['id']][idx2D], appendID = "_2D")  
    fastq_2D <- fq_2D$fastq
  }
  
  ## We update the individual strands to indicate if they are part of a full 2D read
  template <- as_tibble(cbind(template, full_2D = template[['id']] %in% idx2D))
  complement <- as_tibble(cbind(complement, full_2D = complement[['id']] %in% idx2D))
  
  # ## combine the template, complement and 2D data
  baseCalled <- rbind(template, complement)
  fastq <- ShortRead::append(fastq_template, fastq_complement)
  if(length(idx2D)) {
    fastq <- ShortRead::append(fastq, fastq_2D)
  }
  
  message("Done")
  obj <- new("Fast5Summary", 
             readInfo = readInfo, 
             rawData = tibble(), 
             eventData = rawEventData, 
             baseCalled = baseCalled, 
             fastq = fastq,
             versions = list('IONiseR' = strsplit(as.character(packageVersion("IONiseR")),
                                                  ".",
                                                  fixed=TRUE)[[1]]
                             # 'MinKNOW' = max(versions))
             )
  )
  
  return(obj)
}