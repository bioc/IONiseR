## produce a ShorReadQ object from the fast5 string
#' @importFrom Biostrings BStringSet DNAStringSet
#' @importFrom ShortRead FastqQuality
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

## produce a ShorReadQ object from a vector of fast5 strings
#' @importFrom Biostrings BStringSet DNAStringSet
#' @importFrom ShortRead FastqQuality
#' @importFrom BiocGenerics width
.processFastqVec <- function(strings, readIDs = NULL, appendID = NULL) {
    
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


.getFastqString <- function(file, strand = "template") {
  ## returns the unprocess fastq string stored in fast5 files
  fid <- H5Fopen(file)  
  on.exit(H5Fclose(fid))
  
  ## is there any template data?
  exists <- .groupExistsObj(fid, group = paste0("/Analyses/Basecall_2D_000/BaseCalled_", strand))
  if(exists) {
    gid <- H5Gopen(fid, paste0("/Analyses/Basecall_2D_000/BaseCalled_", strand))
    did <- H5Dopen(gid, "Fastq")
    fastq <- H5Dread(did)
    
    H5Dclose(did)
    H5Gclose(gid)
  }
  
  return(fastq)
}
