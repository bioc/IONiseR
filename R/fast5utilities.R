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

## try to guess version of fast5 file
#' @importFrom stringr str_detect
.fast5version <- function(file) {
    
    groups <- unique(h5ls(file)[,'group'])
    
    ## classic version with 2D slot
    if(any(str_detect("/Analyses/Basecall_2D_[0-9]{3}/BaseCalled_template", string = groups))) {
        version <- 1
        ## not often seen variant, seen near r9 chemistry release    
    } else if(any(str_detect("/Analyses/Basecall_RNN_1D_[0-9]{3}/BaseCalled_template", string = groups))) {
        version <- 3
        ## version produced since metrichor 1.16    
    } else if(any(str_detect("/Analyses/Basecall_1D_[0-9]{3}/BaseCalled_template", string = groups))) {
        version <- 2
        ## if none of the above are found, no basecalling has been performed    
    } else {
        version <- 0
    }
    
    return(version)
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