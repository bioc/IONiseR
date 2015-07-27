#' Plot the proportion of template, complement and 2D reads found a dataset.
#' 
#' Generates a bar plot showing the breakdown of read types found in a set of fast5 files.  There is a strict hierarchy to the types of read that can be found in a fast5 file.  A full 2D read requires both a complement and template strand to have been read correctly.  Similarly, a complement strand can only be present if the template was read successfully.  Finally, you can encounter a file containing now called bases on either strand.
#' Here we visualise the total number of fast5 files, along with the counts containing each of the categories above.  For an ideal dataset all four bars will be the same height.  This is unlikely, but the drop between bars can give some indication of data quality.
#' @param summaryData Object of class \linkS4class{Fast5Summary}.
#' @return Returns an object of class \code{gg} representing the plot.
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    plotReadCategoryCounts( s.typhi.rep2 )
#' }
#' @export
#' @importFrom dplyr summarise count
plotReadCategoryCounts <- function(summaryData) {
    tab <- c(nrow(readInfo(summaryData)),
             count(baseCalled(summaryData), strand)[,n],
             count(baseCalled(summaryData), full_2D)[2,n] / 2)
    
    res <- data.table(category = factor(c('Fast5 File Count', 'Template', 'Complement', 'Full 2D'),
                                        levels = c('Fast5 File Count', 'Template', 'Complement', 'Full 2D')),   
                      count = tab)
    
    ggplot(res, aes(x = category, y = count, fill = category)) + 
        geom_bar(stat = 'identity') +
        xlab('read type') +
        guides(fill=FALSE)
}

#' Visualise the mean base quality of each read
#' 
#' Generates a box plot showing the mean base quality for each read, broken down into the three categories of read type that can be found in a fast5 file.
#' 
#' @param summaryData Object of class \linkS4class{Fast5Summary}.
#' @return Returns an object of class \code{gg} representing the plot.
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    plotReadCategoryQuals( s.typhi.rep2 )
#' }
#' @export
#' @importFrom ShortRead alphabetScore
#' @importFrom Biostrings quality
plotReadCategoryQuals <- function(summaryData) {
    
    fq <- fastq(summaryData)
    ## i've inserted an extra level in the factor here as I like the 4 colour option more - 
    ## probaby not the best thing to do.
    readType <- factor(.readtypeFromFASTQ(fq), levels = c('space', 'template', 'complement', '2D'))
    meanBaseQuality <- alphabetScore(quality(fq)) / width(fq) 
    
    res <- data.table(readType, meanBaseQuality)
    ggplot(res, aes(x = readType, y = meanBaseQuality, fill = factor(readType))) + 
        geom_boxplot() +
        xlab('read type') +
        ylab('mean base quality per read') +
        guides(fill=FALSE)
}

#' Plot the number of active channels for each minute of run time
#' @param summaryData Object of class \linkS4class{Fast5Summary}.
#' @param title Character string specifying the plot title.
#' @return Returns an object of class \code{gg} representing the plot.
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    plotActiveChannels( s.typhi.rep2 )
#' }
#' @export
#' @importFrom dplyr summarise count
plotActiveChannels <- function(summaryData, title = "") {
    startEndSummary <- summarise(rawData(summaryData), first = start_time %/% 60, last = (start_time + duration) %/% 60)
    tab <- data.frame(minute = unlist(apply(startEndSummary, 1, function(x) { x[1]:x[2] }))) %>%
        count(minute)
    
    ggplot(tab, aes(x = minute, y = n)) + 
        geom_point() + 
        ylab("Active Channels") +
        ggtitle(title)
}

#' Plot the accumulation of reads over the duration of the experiment.
#' @param summaryData Object of class \linkS4class{Fast5Summary}.
#' @param title Character string specifying the plot title.
#' @return Returns an object of class \code{gg} representing the plot.
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    plotReadAccumulation( s.typhi.rep2 )
#' }
#' @export
#' @importFrom dplyr group_by summarise mutate order_by with_order n
plotReadAccumulation <- function(summaryData, title = "") {
    readAccumulation <- group_by(rawData(summaryData), minute = start_time %/% 60) %>%
        summarise(new_reads = n()) %>%
        mutate(accumulation = order_by(minute, cumsum(new_reads)))
    ggplot(readAccumulation, aes(x = minute, y = accumulation)) + 
        geom_point() + 
        ylab("reads produced") +
        ggtitle(title)
}


#' Plot the accumulation of reads over the duration of the experiment.
#' @param summaryData Object of class \linkS4class{Fast5Summary}.
#' @param title Character string specifying the plot title.
#' @return Returns an object of class \code{gg} representing the plot.
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    plotReadAccumulation( s.typhi.rep2 )
#' }
#' @export
#' @importFrom dplyr group_by summarise mutate order_by with_order n
plotYield <- function(summaryData, title = "") {
    readAccumulation <- group_by(rawData(summaryData), minute = start_time %/% 60) %>%
        summarise(new_reads = n()) %>%
        mutate(accumulation = order_by(minute, cumsum(new_reads)))
    ggplot(readAccumulation, aes(x = minute, y = accumulation)) + 
        geom_point() + 
        ylab("reads produced") +
        ggtitle(title)
}

#' Plot the mean rate at which events occur 
#' 
#' For each read, the ratio between the number of events comprising the read and the time spent in the pore is calculated.  This is then plotted against the time the read entered the pore, allow us to assess whether the rate at which events occur changes during the experiment run time.
#' @param summaryData Object of class \linkS4class{Fast5Summary}.
#' @return Returns an object of class \code{gg} representing the plot.
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    plotEventRate( s.typhi.rep2 )
#' }
#' @export
plotEventRate <- function(summaryData) {
    ggplot(rawData(summaryData), aes(x = start_time %/% 60, y = num_events / duration)) + 
        geom_point(alpha = 0.3) + 
        ylim(0,80) +
        xlab("minute")
}

#' Plot the mean rate at which bases are recorded 
#' 
#' For each read, the ratio between the total number of bases called in the read
#' (template and complement strand, but not 2D composite) and the time spent in
#' the pore is calculated.  This is then plotted against the time the read 
#' entered the pore, allow us to assess whether the rate at which callable 
#' bases are read changes during the experiment run time.
#' 
#' This is likely very similar to \code{\link{plotEventRate}}, although one may
#' find that large number of events occur that can not be base called, 
#' resulting in a difference between these two plots.
#' 
#' 
#' @param summaryData Object of class \linkS4class{Fast5Summary}.
#' @return Returns an object of class \code{gg} representing the plot.
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    plotBaseProductionRate( s.typhi.rep2 )
#' }
#' @export
#' @importFrom dplyr mutate 
#' @importFrom ShortRead width
plotBaseProductionRate <- function(summaryData) {
    ## select only the fastq records for the template and complement reads
    ## ignore the composite 2D reads here
    fastqIDX <- .matchRecords(summaryData)[,c(fastqTemplate, fastqComplement)]
    fastqIDX <- fastqIDX[-which(is.na(fastqIDX))]
    res <- mutate(baseCalled(summaryData), bases_called = width(summaryData@fastq[ fastqIDX ]))
    ggplot(res, aes(x = start_time %/% 60, y = bases_called / duration)) + 
        geom_point(alpha = 0.3) + 
        ylim(0,80) +
        xlab("minute")
}

#' View changes in signal against run time.
#' 
#' Plots the median recorded current for each fast5 file against the time at 
#' which the recording began.
#' 
#' @param summaryData Object of class \linkS4class{Fast5Summary}.
#' @return Returns an object of class \code{gg} representing the plot.
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    plotCurrentByTime( s.typhi.rep2 )
#' }
#' @export
plotCurrentByTime <- function(summaryData) {
    ggplot(rawData(summaryData), aes(x = start_time, y = median_signal)) + 
        geom_point() + 
        xlab("time (seconds)") +
        ylab("current (pA)")
}



readTypesByTime <- function(summaryData, minute_group = 10) {
    tmp <- left_join(baseCalled(summaryData), readInfo(summaryData), by = 'id') %>%
        filter(strand == "template") %>%
        group_by(time_group = start_time %/% (60 * minute_group), full_2D, pass) %>%
        summarise(count = n(), hour = (time_group * minute_group)/60 )
    
    ggplot(tmp, aes(x = hour, y = count, colour = interaction(full_2D, pass))) + 
        geom_point(size = 3) +
    scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9"), 
                      name="Category",
                      labels=c("Not 2D", "2D - Fail", "2D - Pass"))
}


