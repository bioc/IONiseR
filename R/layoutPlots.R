#' Create layout plot of flowcell
#' 
#' Creates a plot representing the layout of a MinION flow cell.  Each circle
#' represents an individual channel with the intensity relecting a specified 
#' sequencing metric.  This function is a more generalised version of 
#' \code{\link{layoutPlot}}, allowing the user to map any value the like on 
#' the channel layout.
#' @param data A data.frame or data.table.  Should have at least two columns, 
#' one of which has the name 'channel'.
#' @param zValue Character string specifying the name of the column to be used 
#' for the colour scaling.
#' @return Returns an object of \code{gg} representing the plot.
#' @examples
#' library(dplyr)
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    ## calculate and plot the mean number of events recorded by each channel
#'    avgEvents <- left_join(readInfo(s.typhi.rep2), rawData(s.typhi.rep2), by = 'id') %>% 
#'    group_by(channel) %>% 
#'    summarise(mean_nevents = mean(num_events))
#'    channelHeatmap(avgEvents, zValue = 'mean_nevents')
#' }
#' @export
channelHeatmap <- function(data, zValue) {
    
    if(!zValue %in% colnames(data)) {
        stop("Column '", zValue, "' not found")
    }
    
    plottingMap <- data.table(.channelToXY(data[,channel], shiftRows = TRUE), zValue = data[,get(zValue)])
    
    ggplot(plottingMap, aes_string(x = "matrixCol", y = "matrixRow", colour = "zValue")) + 
        geom_rect(mapping = aes(xmin = 0, xmax = 8.5, ymin = -0.5, ymax = 33.5), fill = "grey50", colour = "black") +
        geom_rect(mapping = aes(xmin = 9, xmax = 17.5, ymin = -0.5, ymax = 33.5), fill = "grey50", colour = "black") +
        geom_point(size = 5.5, colour = "black") +
        geom_point(size = 4.5) + 
        scale_colour_gradient(low="darkblue", high="green") + 
        theme(panel.background = element_rect(fill = "grey50"), 
              panel.grid.major = element_line(colour = "grey50"), 
              panel.grid.minor = element_line(colour = "grey50")) + 
        scale_x_continuous(breaks=NULL) +
        scale_y_continuous(breaks=NULL) +
        xlab("") + 
        ylab("")
}

#' Create layout plot of flowcell
#' 
#' Creates a plot representing the layout of a MinION flow cell.  Each circle represents an idividual channel with the intensity relecting the total kilobases of sequence produced.  This only considers reads marked as template or complement, 2D reads are ignored as they are generated from the former two.
#' @param summaryData Object of class \linkS4class{Fast5Summary}.
#' @param attribute Character string indicating what to plot. Currently accepted values are: "nreads", "kb", "signal".
#' @return Returns an object of \code{gg} representing the plot.
#' @examples
#' if( require(minionSummaryData) ) {
#'    data(s.typhi.rep2, package = 'minionSummaryData')
#'    layoutPlot( s.typhi.rep2, attribute = 'nreads' )
#'    layoutPlot( s.typhi.rep2, attribute = 'kb' )
#' }
#' @export
#' @importFrom dplyr group_by summarise n
#' @importFrom BiocGenerics width
layoutPlot <- function(summaryData, attribute = NULL) {
    
    if(attribute == "nreads") {
        res <- group_by(summaryData@readInfo, channel) %>% summarise(nreads = n())
    } else if (attribute == "kb") {
        res <- mutate(baseCalled(summaryData), 
                      seq_length = width(summaryData@fastq[1:nrow(summaryData@baseCalled)]), 
                      channel = readInfo(summaryData)[ match(baseCalled(summaryData)[,id], readInfo(summaryData)[,id]), channel]) %>%
            group_by(channel) %>%
            summarise(kb = sum(seq_length) / 1000)  
    } else if (attribute == "signal") {
        res <- mutate(summaryData@rawData, channel = summaryData@readInfo[match(summaryData@rawData[,id], summaryData@readInfo[,id]), channel]) %>% 
            group_by(channel) %>% 
            summarise(signal = mean(median_signal))
    } else {
        error("Unsupported attribute")
    }
    
    channelHeatmap(res, zValue = attribute) 
}


.channelToXY <- function(pore = 1:512, insertGap = TRUE, shiftRows = TRUE) {
    
    pore <- pore-1
    
    side <- (pore %/% 32) %% 2
    
    rowblock <- pore %/% 64
    rowblock[rowblock > 0] <- abs(rowblock[rowblock > 0] - 8)
    
    row <- (rowblock * 4) + ((pore %% 32) %/% 8) + 1
    column <- (-(side-1) * (pore %% 8)) + (side  * (15 - pore %% 8)) + 1
    
    if(insertGap) {
        column[column > 8] <- column[column > 8] + 1
    }
    
    if(shiftRows) {
        column[ (which(!row %% 2)) ] <- column[ (which(!row %% 2)) ] - 0.5
    }
    
    return(data.frame(channel = pore+1, matrixRow = row, matrixCol = column))
}


