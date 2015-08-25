#' Display correlation between pentemer proportions in two time windows
#' 
#' Plots 
#' 
#' @param summaryData Object of class \linkS4class{Fast5Summary}.
#' @param kmerLength Specifies the length of kmers to compare. Defaults to 5 
#' given the current pentamer reading nature of the nanopores.
#' @param groupedMinutes Defines how many minutes each grouping of reads spans.
#' @param only2D  Logical. If TRUE kmers are computed for only full 2D reads.
#' If FALSE 2D reads are ignored and all available template and complement 
#' strands are used.
#' @return Returns an object of class \code{gg} representing the plot.
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep3, package = 'minionSummaryData')
#'    plotKmerFrequencyCorrelation( s.typhi.rep3, only2D = FALSE )
#' }
#' @export
#' @importFrom dplyr summarise_each funs arrange
#' @importFrom Biostrings oligonucleotideFrequency
#' @importFrom tidyr spread gather
#' @importFrom ShortRead sread
plotKmerFrequencyCorrelation <- function(summaryData, kmerLength = 5, groupedMinutes = 10, only2D = TRUE) {
    
    if(only2D) {
        idx <- grep('2D', id(fastq(summaryData)))
    } else {
        idx <- 1:nrow(baseCalled(summaryData))
    }
    pentamers <- oligonucleotideFrequency(x = sread(fastq(summaryData)[idx]), width = kmerLength, as.prob = TRUE)
    
    if(only2D) {
        tmp <- data.table(filter(baseCalled(summaryData), strand == "template", full_2D == TRUE), pentamers)
    } else {
        tmp <- data.table(baseCalled(summaryData), pentamers)
    }
    
    tmp2 <- group_by(tmp, time_group = start_time %/% (60 * groupedMinutes)) %>%
        summarise_each(funs(mean), AAAAA:TTTTT) %>%
        arrange(time_group) %>%
        gather(key = "pentamer", value = "freq", AAAAA:TTTTT) %>%
        spread(key = time_group, value = freq) %>% 
        select(-pentamer)
    
    correlations <- cor(as.matrix(tmp2))
    
    hours <- (as.integer(colnames(correlations)) * groupedMinutes) / 60
    tmp3 <- data.table(x = rep(hours, each = length(hours)), y = rep(hours, length(hours)), cor = as.vector(correlations))
    
    ggplot(tmp3, aes(x = x, y = y, fill = cor)) + 
        geom_raster() +
        scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) +
        xlab("hour") +
        ylab("hour")
}
