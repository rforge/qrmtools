### Graphics for displaying data ###############################################

##' @title Image displaying NAs in a data set
##' @param x A matrix (ideally an xts object)
##' @param col The colors for NA and non-NA, respectively
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param text See mtext()
##' @param side See mtext()
##' @param line See mtext()
##' @param adj See mtext()
##' @param ... Additional arguments passed to image()
##' @return invisible()
##' @author Marius Hofert
plot_NA <- function(x, col = c("black", "white"), xlab = "Time", ylab = "Component",
                    text = "Black: NA; White: Available data", side = 4, line = 1, adj = 0,
                    ...)
{
    stopifnot(is.matrix(x))
    x. <- if(inherits(x, "xts")) {
        index(x) # use the time points
    } else {
        rn <- rownames(x)
        if(is.null(rn)) seq_len(nrow(x)) else rn # if available, use row names, otherwise numbers
    }
    image(x = x., y = seq_len(ncol(x)), z = is.na(x),
          col = rev(col), xlab = xlab, ylab = ylab, ...)
    if(!is.null(text) && nchar(text) > 0)
        mtext(text, side = side, line = line, adj = adj)
    invisible()
}

##' @title Plot of a Matrix
##' @param x A matrix
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param scales See levelplot(); if NULL, labels and ticks are omitted
##' @param at See levelplot()
##' @param colorkey See levelplot()
##' @param col.regions See levelplot()
##' @param ... Additional arguments passed to levelplot()
##' @return The level plot
##' @author Marius Hofert
##' @note Another option would be:
##'       corrplot::corrplot(err, method="color", col=grey(seq(0.4, 1, length.out=200)),
##'                          tl.col="black", is.corr=FALSE)
plot_matrix <- function(x, xlab="Column", ylab="Row",
                        scales=list(alternating=c(1,1), tck=c(1,0), x=list(rot=90)),
                        at=seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=200),
                        colorkey=list(at=at),
                        col.regions=colorRampPalette(c("royalblue3", "white", "maroon3"),
                                                     space="Lab")(200),
                        ...)
{
    stopifnot(is.matrix(x), (nr <- nrow(x)) >= 1,
              is.null(scales) || is.list(scales), is.numeric(at),
              is.list(colorkey))
    ran <- range(x, na.rm=TRUE)
    if(min(x, na.rm=TRUE) < at[1] || max(x, na.rm=TRUE) > at[length(at)])
        stop("'x' values outside the range spanned by 'at'. Choose 'at' appropriately.")
    levelplot(t(x)[,nr:1], xlab=xlab, ylab=ylab,
              col.regions=col.regions,
              scales=if(is.null(scales)) list(alternating=c(0,0), tck=c(0,0)) else scales,
              at=at, colorkey=colorkey, ...)
}
