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
    
    if("pass" %in% names(readInfo(summaryData))) {
    pf <- any(count(readInfo(summaryData), pass)[,pass], na.rm = TRUE)
    } else {
        pf <- FALSE
    }
    if(pf) {
        tab <- c(tab, pass[pass==TRUE,n])
        res <- data.table(
            category = factor(c('Fast5 File Count', 'Template', 'Complement', 'Full 2D', 'Pass'),
                              levels = c('Fast5 File Count', 'Template', 'Complement', 'Full 2D', 'Pass')),   
            count = tab)
    } else {
        res <- data.table(
            category = factor(c('Fast5 File Count', 'Template', 'Complement', 'Full 2D'),
                              levels = c('Fast5 File Count', 'Template', 'Complement', 'Full 2D')),   
            count = tab)
    }

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
#' @return Returns an object of class \code{gg} representing the plot.
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    plotActiveChannels( s.typhi.rep2 )
#' }
#' @export
#' @importFrom dplyr summarise count
plotActiveChannels <- function(summaryData) {
    startEndSummary <- summarise(rawData(summaryData), first = start_time %/% 60, last = (start_time + duration) %/% 60)
    tab <- data.frame(minute = unlist(apply(startEndSummary, 1, function(x) { x[1]:x[2] }))) %>%
        count(minute)
    
    ggplot(tab, aes(x = minute / 60, y = n)) + 
        geom_point() + 
        xlab("hour") +
        ylab("active channels")
}

#' Plot the accumulation of reads over the duration of the experiment.
#' @param summaryData Object of class \linkS4class{Fast5Summary}.
#' @return Returns an object of class \code{gg} representing the plot.
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    plotReadAccumulation( s.typhi.rep2 )
#' }
#' @export
#' @importFrom dplyr group_by summarise mutate order_by with_order n
plotReadAccumulation <- function(summaryData) {
    readAccumulation <- group_by(rawData(summaryData), minute = start_time %/% 60) %>%
        summarise(new_reads = n()) %>%
        mutate(accumulation = order_by(minute, cumsum(new_reads)))
    ggplot(readAccumulation, aes(x = minute / 60, y = accumulation)) + 
        geom_point() + 
        xlab("hour") +
        ylab("reads produced")
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
    ggplot(rawData(summaryData), aes(x = start_time / 3600, y = median_signal)) + 
        geom_point() + 
        xlab("hour") +
        ylab("median current (pA)")
}



plotReadTypesByTime <- function(summaryData, groupedMinutes = 10) {
    tmp <- left_join(baseCalled(summaryData), readInfo(summaryData), by = 'id') %>%
        filter(strand == "template") %>%
        group_by(time_group = start_time %/% (60 * groupedMinutes), full_2D, pass) %>%
        summarise(count = n(), hour = (time_group * groupedMinutes)/60 )
    
    ggplot(tmp, aes(x = hour, y = count, colour = interaction(full_2D, pass))) + 
        geom_point(size = 3) +
    scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9"), 
                      name="Category",
                      labels=c("Not 2D", "2D - Fail", "2D - Pass"))
}


plot2DYield <- function(summaryData, groupedMinutes = 1) {
    
    only2d <- .get2D(summaryData)
    tmp.fq <- fastq(only2d)[grep("2D", id(fastq(only2d))),]
    
    tmp <- inner_join(readInfo(only2d), data.table(id = IONiseR:::.idFromFASTQ(tmp.fq), nbases = width(tmp.fq)), by = "id")
    tmp <- inner_join(tmp, baseCalled(only2d), by = "id")
    
    readAccumulation <- group_by(tmp, time_group = start_time %/% (60 * groupedMinutes), pass) %>%
        summarise(nbases = sum(nbases), hour = (time_group * groupedMinutes)/60 ) %>%
        arrange(time_group) %>%
        group_by(pass) %>%
        mutate(accumulation = order_by(time_group, cumsum(nbases)))
    ggplot(readAccumulation, aes(x = hour, y = accumulation, colour = pass)) + 
        geom_point() + 
        ylab("bases produced")
}


channelActivityPlot <- function(summaryData, zScale = NULL, zAverage = TRUE) {
    
    tmp <- left_join(readInfo(summaryData), rawData(summaryData), by = 'id') %>%
        select(id, channel, start_time, duration)
    
    ## if we've provided a zvalue, add it to our table and set the column name
    if(!is.null(zScale)) {
        setnames(zScale, names(zScale)[ncol(zScale)], "zvalue")
        tmp <- left_join(tmp, zScale, by = "id")
    } else {
        tmp <- mutate(tmp, zvalue = 1)
    }
    
    p1 <- ggplot(NULL) +
        geom_segment(data = tmp,
                     mapping = aes(y = channel, yend = channel, 
                                   x = (start_time / 3600), xend = (start_time + duration) / 3600,
                                   color = zvalue), 
                     size = 1) +
        xlab("hour") +
        scale_x_continuous(expand = c(0, 0))
    
    if(zAverage) {
        tmp2 <- tmp %>% 
            group_by(time_bin = start_time %/% 600) %>%
            summarise(mean_value = mean(zvalue))
        p1 <- p1 + geom_rect(mutate(tmp2, start_time = (time_bin * 600) / 3600), mapping = aes(xmin = start_time, xmax = start_time + (600/3600), ymin = -100, ymax = -25, fill = mean_value))
    }
        
    p1
}
