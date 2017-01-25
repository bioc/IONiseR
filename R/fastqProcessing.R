## produce a ShorReadQ object from the fast5 string
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
#' @importFrom stringr str_length
.processFastqVec <- function(strings, readIDs = NULL, appendID = NULL) {
    
    if(any(str_length(strings) == 0))
        strings <- strings[-which(str_length(strings) == 0)]
    
    strings <- str_replace(string = strings, pattern = "^@", replacement = "")
    
    fastqStrings <- strsplit(strings, "\n")
    fastqStrings <- do.call(rbind, fastqStrings)
    
    if(is.null(readIDs)) {
        id <- BStringSet(paste0(fastqStrings[,1], appendID))
    } else {
        id <- BStringSet(paste0(readIDs, appendID))
    }
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
  ## returns the unprocess fastq string stored in fast5 files
    if(is.character(file)) {
        fid <- H5Fopen(file)
        on.exit(H5Fclose(fid))
    } else {
        fid <- file
    }

    if( dontCheck || .groupExistsObj(fid, group = paste0("/Analyses/Basecall_", d, "_000/BaseCalled_", strand, "/Fastq")) ) {
        #gid <- H5Gopen(fid, paste0("/Analyses/Basecall_", d, "_000/BaseCalled_", strand))
        did <- H5Dopen(fid, paste0("/Analyses/Basecall_", d, "_000/BaseCalled_", strand, "/Fastq"))
       
        fastq <- H5Dread(did)
        
        H5Dclose(did)
        #H5Gclose(gid)
    } else {
        fastq <- ""
    }
    
    return(fastq)
}

.processStrandSpecifier <- function(strand) {
    
    strand <- stringr::str_split(strand, pattern = "\\|", simplify = TRUE)[1,]
    if(any(!strand %in% c("template", "complement", "2D", "all", "both"))) {
        stop("Unexpected option provided to 'strand'")
    }
    
    if("all" %in% strand) {
        strand <- c("template", "complement", "2D")
    } else if ("both" %in% strand) {
        strand <- c("template", "complement")
    }
    
    return(strand)
}


#' @export
#' @importFrom ShortRead writeFastq
fast5toFastq <- function(files, strand = "all", fileName = NULL, 
                         outputDir = NULL, yield = NULL, dontCheck = TRUE) {
    
    ## understand the file structure
    status <- IONiseR:::.fast5status(files = sample(files, size = min(length(files), 15)))
    
    strand <- .processStrandSpecifier(strand)
    ## if only 1D pipeline has been run, we can't get complement/2d data
    if(!status$basecall_2d && any(c("complement", "2D") %in% strand)) {
        warning("This data has only been processed using the 1D pipeline.\n",
                "Only FASTQ files for the template strand will be generated",
                call. = FALSE)
        strand <- "template"
    }
    
    for(s in strand) {
        fastq_vec <- bpmapply(FUN = IONiseR:::.getFastqString, files, 
                            strand = s, dontCheck = dontCheck,
                            USE.NAMES = FALSE)
        fastq_fq <- IONiseR:::.processFastqVec(fastq_vec)$fastq
        ShortRead::writeFastq(fastq_fq, 
                              file = file.path(outputDir, paste0(fileName, "_", s, ".fq.gz")))
    }
}
