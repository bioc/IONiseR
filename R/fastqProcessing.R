## produce a ShortReadQ object from the fast5 string
#' @importFrom Biostrings BStringSet DNAStringSet
#' @importFrom ShortRead FastqQuality
#' @importFrom methods new
.processFastq <- function(string, readID = NULL, appendID = NULL) {
    
    fastqString <- strsplit(string, "\n")[[1]]
    if(is.null(readID)) {
        id <- BStringSet(paste0(fastqString[1], appendID))
    } else {
        id <- BStringSet(paste0(readIDs, appendID))
    }
    sread <- DNAStringSet(fastqString[2])
    qual <- FastqQuality(fastqString[4])
    
    new("ShortReadQ", id = id, sread = sread, quality = qual)
}

## produce a ShortReadQ object from a vector of fast5 strings
#' @importFrom Biostrings BStringSet DNAStringSet
#' @importFrom ShortRead FastqQuality
#' @importFrom BiocGenerics width
#' @importFrom stringr str_length str_match str_replace_all
.processFastqVec <- function(strings, readIDs = NULL, appendID = NULL) {
    
    if(any(str_length(strings) == 0))
        strings <- strings[-which(str_length(strings) == 0)]
    
    strings <- str_replace(string = strings, pattern = "^@", replacement = "")
    
    #fastqStrings <- strsplit(strings, "\n")
    #fastqStrings <- do.call(rbind, fastqStrings)
    fastqStrings <- stringr::str_match(strings, pattern = "(.*)\n(.*)\n(\\+)\n(.*)")[,2:5,drop=FALSE]
    
    if(is.null(readIDs)) {
        id <- BStringSet(paste0(fastqStrings[,1], appendID))
    } else {
        id <- BStringSet(paste0(readIDs, appendID))
    }
    
    ## replace U with T to cope with RNA data
    ## this could be made more generalisable
    fastqStrings[,2] <- stringr::str_replace_all(string = fastqStrings[,2], 
                             pattern = "U", 
                             replacement = "T")
    
    sread <- DNAStringSet(fastqStrings[,2])
    qual <- FastqQuality(fastqStrings[,4])
    
    ## occasionally we encounter records where 
    ## the number of bases != the number of qualities
    ## these are not valid and should be removed
    invalid <- which(width(sread) != width(qual))
    if(length(invalid)) {
        id <- id[-invalid]
        sread <- sread[-invalid]
        qual <- qual[-invalid]
    }
    
    fastq <- new("ShortReadQ", id = id, sread = sread, quality = qual)
    return(list(fastq = fastq, invalid = invalid))
}

#' @importFrom stringr str_replace
.getFastqString <- function(file, strand = "template", d = "1D",
                            dontCheck = TRUE) {
  ## returns the unprocessed fastq string stored in fast5 files
    if(is.character(file)) {
        fid <- H5Fopen(file)
        on.exit(H5Fclose(fid))
    } else {
        fid <- file
    }

    group <- paste0("/Analyses/Basecall_", d, "_000/BaseCalled_", strand, "/Fastq")
    if( dontCheck || .groupExistsObj(fid, group = group) ) {
        did <- H5Dopen(fid, group)
        fastq <- H5Dread(did)
        H5Dclose(did)
    } else {
        fastq <- ""
    }
    
    return(fastq)
}

.processStrandSpecifier <- function(strand) {
    
    strand <- stringr::str_split(strand, pattern = "\\|", simplify = TRUE)[1,]
    if(any(!strand %in% c("template", "complement", "2D", "all", "both"))) {
        stop("Unexpected option provided to 'strand'.",
             "Please provide some combination of the following options:",
             "'template', 'complement', '2D', 'all', 'both'")
    }
    
    if("all" %in% strand) {
        strand <- c("template", "complement", "2D")
    } else if ("both" %in% strand) {
        strand <- c("template", "complement")
    }
    
    return(strand)
}


#' Extract FASTQ files from fast5 files
#' 
#' This function provides direct access to the FASTQ entries held within fast5
#' files.  If you are only interested in getting hold of the base called reads,
#' and don't require any raw-signal or event information, use this function.
#' Given a vector of fast5 files, the FASTQ entries will be combined and up to
#' three gzip compressed FASTQ will be created - one for each of the template, 
#' complement and 2D strands depending upon what is available in the input 
#' files.
#' 
#' @param files Character vector of fast5 files to be read.
#' @param strand Character vector specifying the strand to extract.  Can take
#' any combination of the following options: "template", "complement", "2D", 
#' "all", "both".
#' @param fileName Stem for the name of the names of the output file names.
#' The appropriate strand will be appended to each file e.g. 
#' fileName_complement.fq.gz or fileName_template.fq.gz
#' @param outputDir Directory output files should be written to.
#' @param ncores Specify the number of CPU cores that should be used to process 
#' the files.  Currently this seems to be more IO bound than CPU, so there
#' is little benefit achieved by using a high number of cores.
#' 
#' @return No value returned.  Run for the side effect of writing the FASTQ
#' files to disk.
#' 
#' @examples \dontrun{
#' fast5files <- list.files('/foo/bar/', pattern = '.fast5$')
#' summaryData <- readFast5Summary(fast5files)
#' }
#' 
#' @export
#' @importFrom ShortRead writeFastq
#' @importFrom BiocParallel MulticoreParam bpmapply register
fast5toFastq <- function(files, strand = "all", fileName = NULL, 
                         outputDir = NULL, ncores = 1) {
    
    ## check files exist and can be accessed
    files <- files[which(file.exists(files))]
    if(length(files) == 0) {
        stop('None of the provided files can be accessed. ',
            'Have you supplied the correct path?')
    }
    
    ## understand the file structure
    status <- .fast5status(files = sample(files, size = min(length(files), 15)))
    
    strand <- .processStrandSpecifier(strand)
    ## if only 1D pipeline has been run, we can't get complement/2d data
        if( (!nchar(status$loc_complement) && "complement" %in% strand) ) {
        warning("This data has only been processed using the 1D pipeline.\n",
                "FASTQ files for the complement strand will not be generated.",
                call. = FALSE)
            strand <- strand[-which(strand == "complement")]
        }
    
    ## if only 1D pipeline has been run, we can't get complement/2d data
    if( !status$basecalled_2d && "2D" %in% strand ) {
        warning("This data does not contain 2D base calls.\n",
                "FASTQ files for the 2D strand will not be generated.",
                call. = FALSE)
        strand <- strand[-which(strand == "2D")]
    }
    
    register(MulticoreParam(workers = ncores))
    
    for(s in strand) {
        d <- str_match(status[[ paste0('loc_', s) ]], 
                       pattern = "Basecall_([12]D)")[,2]
        fastq_vec <- bpmapply(FUN = .getFastqString, files, 
                            MoreArgs = list(strand = s, d = d, dontCheck = FALSE),
                            USE.NAMES = FALSE)
        if(length(fastq_vec) == 0) {
            warning("No entries found for ", s, "strand. ",
                    "No FASTQ file created.",
                    call. = FALSE)
            next()
        }
        fastq_fq <- .processFastqVec(fastq_vec)$fastq
        writeFastq(fastq_fq, 
                   file = file.path(outputDir, 
                                    paste0(fileName, "_", s, ".fq.gz")))
    }
}
