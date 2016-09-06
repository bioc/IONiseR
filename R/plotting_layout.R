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
#' @importFrom dplyr as_data_frame
#' @export
channelHeatmap <- function(data, zValue) {
    
    if(!zValue %in% colnames(data)) {
        stop("Column '", zValue, "' not found")
    }
    
    plottingMap <- as_data_frame(cbind(.channelToXY(data[['channel']], shiftRows = TRUE), zValue = data[[zValue]]))
    
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
#' @export
muxHeatmap <- function(data, zValue) {
    
    if(!zValue %in% colnames(data)) {
        stop("Column '", zValue, "' not found")
    }
    
    tmp <- .muxToXY(.channelToXY(shiftRows = FALSE, insertGap = FALSE), shiftRows = TRUE)
    plottingMap <- data.table(right_join(tmp, data, by = c("mux" = "mux", "channel" = "channel")), zValue = data[[zValue]])
    
    tmp <- left_join(tmp, group_by(plottingMap, channel) %>% summarise(zValue = sum(zValue)), by = "channel")

    channel_data <- group_by(tmp, channel) %>% 
        summarise(row = mean(matrixRow), col = mean(matrixCol), meanZValue = mean(zValue))
    
    ggplot() +
       # geom_rect(mapping = aes(xmin = -0.5, xmax = 33, ymin = -0.5, ymax = 33.5), fill = "grey50", colour = "black") +
      #  geom_rect(mapping = aes(xmin = 35.5, xmax = 69, ymin = -0.5, ymax = 33.5), fill = "grey50", colour = "black") +
        geom_rect(data = channel_data, 
                  mapping = aes(xmax = col+2, xmin = col-2, ymax = row+0.5, ymin = row-0.5, fill = meanZValue), color = "gray40", alpha = 1) +
        #geom_point(data = plottingMap, mapping = aes(x = matrixCol, y = matrixRow, col = zValue), size = 5) +
        geom_polygon(data = data.table(id = rep(1:nrow(plottingMap), each = 50), 
                                       zValue = rep(plottingMap[,zValue], each = 50),
                                       rbindlist(apply(as.matrix(plottingMap[,4:3, with = FALSE]), 1, circleFun, npoints = 50, diameter = 0.8))), aes(x = x, y = y, group = id, fill = zValue), colour = "grey20") + 
        #scale_colour_gradient(low="darkblue", high="green") + 
        scale_fill_gradient(na.value = "white", low = "#999999", high = "orange") + 
        theme(panel.background = element_rect(fill = "white"), 
              panel.grid.major = element_line(colour = "white"), 
              panel.grid.minor = element_line(colour = "white")) + 
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
        res <- group_by(readInfo(summaryData), channel) %>% summarise(nreads = n())
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

.muxToXY <- function(channelMap, insertGap = TRUE, shiftRows = TRUE) {
    
    muxMap <- data.table(mux = rep(1:4, each = nrow(channelMap)),
                         rbind(channelMap, channelMap, channelMap, channelMap)) %>%
        mutate(oddEven = matrixCol %% 2)

    ##odd columns
    muxMap[mux == 1 & oddEven == 1,matrixCol:=(matrixCol*4)-1]
    muxMap[mux == 2 & oddEven == 1,matrixCol:=(matrixCol*4)]
    muxMap[mux == 3 & oddEven == 1,matrixCol:=(matrixCol*4)-3]
    muxMap[mux == 4 & oddEven == 1,matrixCol:=(matrixCol*4)-2]
    
    ##even columns
    muxMap[mux == 1 & oddEven == 0,matrixCol:=(matrixCol*4)-2]
    muxMap[mux == 2 & oddEven == 0,matrixCol:=(matrixCol*4)-3]
    muxMap[mux == 3 & oddEven == 0,matrixCol:=(matrixCol*4)]
    muxMap[mux == 4 & oddEven == 0,matrixCol:=(matrixCol*4)-1]
    
    if(insertGap) {
        muxMap[matrixCol > 32,matrixCol:=matrixCol+4]
    }
    
    if(shiftRows) {
        muxMap[as.logical(matrixRow %% 2),matrixCol:=matrixCol-0.5]
    }
    
    muxMap <- select(muxMap, -oddEven)
    muxMap
}

.circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}
