#'
#' @title Find Cutoff
#' Finds cells which are positive inputted signal, used in the context of tetramer based CITEseq
#' Make mode tree to calculate peaks
#' Find best bandwidth and get peak values for first two peaks, assuming data is bimodal
#' Find minimum between peaks
#'
#' @param data Dataframe with 1 column containing dex calls and rownames as cell barcodes
#' @param n_peaks Number of peaks expected in data, Default 2 (Bimodal)
#' @param display Prints plot of cutoff
#' @param approach Which approach to use forcalculating threshold
#' @param threhold Manual numericl value, will bypass other approaches
#'
#' @export
#'
#' @import multimode
#' @import mclust
#' @import ggplot2
#' @import tidyverse
#'
#' @return list of two values, first value is a filtered dataframe of values above thresfold and second is ggplot which cutoff used
find_cutoff <- function(data, n_peaks=2, title='',display=TRUE, approach = 1, threshold=NULL){
    #Remove 0's
    data <- data[data[,1] > 0,,drop=FALSE]
    values <- data[,1]
    #Build mode tree
    tree <- multimode::modetree(values, display = FALSE)
    # find lowest band width to contain only n peaks
    best <- tree$locations[rowSums(!is.na(tree$locations)) == n_peaks,]
    # location row.names are lost if filtered data has less than 1 row...
    # different conditions based on datatype of best
    # TODO I think converting locations to dataframe will make if/else statement redundant?
    if("matrix" %in% class(best)){
        best <- best %>% tail(1)

        # get bw and x-intercepts for peaks
        band_width <- row.names(best) %>%
            as.numeric()

        peaks <- best %>%
            as.vector() %>%
            na.omit() %>% sort()
    }

    else{
        row_id <- which(rowSums(!is.na(tree$locations)) == n_peaks)
        band_width <- row.names(tree$locations)[row_id] %>%
            as.numeric()
        peaks <- tree$locations[row_id,] %>%
            as.vector() %>%
            na.omit() %>% sort()
    }

    d <- density(values, bw=band_width)
    d <- data.frame(x=d$x,density=d$y)

    # Best minima
    if(approach == 1){
        # Approach one, find local minimum between two peaks
        min <- d[which(d$x > peaks[1] & d$x < peaks[2]),]
        min <- min[which.min(min$density),'x']
    }
    # Nearest Minima
    if(approach== 2){
        # Approach two, find first local minimum
        min <- d[which(d$x > peaks[1])+1,]
        min <- min[find_local_min(min$density),'x']
    }
    # First derivative
    if(approach==3){
        # Approach three, calculate first derivative,
        # find point which values draw close to 0
        min <- d[which(d$x > peaks[1])+1,]
        slopes <- c()
        for(i in seq_along(min$x)){
            p1 <- c(d$x[i], d$density[i])
            p2 <- c(d$x[i+1], d$density[i+1])

            slope <- (p2[2] - p1[2]) / (p2[1] - p1[1])
            slopes <- append(slopes, atan(slope) * 180 / pi)
        }
        min <- min[which(slopes < -1),'x'][1]
    }
    # Mixed gaussian mode, find values within first peak and calculate
    # mean and 95% CI
    if(approach==4){
        model <- mclust::Mclust(d$x, G=n_peaks, model='V')
        peak1 <- d
        peak1$classification = model$classification
        peak1 <- peak1[peak1$classification == 1,]
        stats <- t.test(peak1$x)
        mean <- stats$estimate
        ci <- stats$conf.int[[2]]
        min <- mean + ci
    }
    # Hard code
    if(approach ==5){min = threshold}
    plot <- ggplot(data, aes(x=values)) +
        geom_density(bw=band_width, linewidth=1) +
        geom_jitter(aes(x=values, y=max(d$density)/2), size=0.5, colour='grey', alpha=0.5) +
        geom_vline(xintercept = peaks, color = "red", size=0.5) +
        geom_vline(xintercept = min, color = "blue", size=0.5) +
        annotate("rect", xmin = min, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.4, fill='sky blue') +
        theme_bw() + ylab("Density") + xlab("Signal") + ggtitle(title)+ xlim(-.5,6)
    if(display){
        print(plot)
    }
    return(list(data[which(values > min),, drop = FALSE], plot))
}


#' @title Local Minima
#' Find location of first local minima
#' within a sequence of numerical values.
find_local_min <- function(sequence){
    for( i in seq_along(sequence)){
        if(sequence[i] < sequence[i+1])
        {
            return(i)
        }
    }
}
