jMincApply <- function (filenames, function.string, mask = NULL, maskval = NULL, reduce = FALSE) 
{
    x <- mincGetVoxel(filenames, 0, 0, 0)
    if (is.null(maskval)) {
        minmask = 1
        maxmask = 99999999
    }
    else {
        minmask = maskval
        maxmask = maskval
    }
    test <- eval(function.string)
    results <- .Call("minc2_model", as.character(filenames), 
        function.string, NULL, as.double(!is.null(mask)), as.character(mask), as.double(minmask), as.double(maxmask), .GlobalEnv, as.double(length(test)), as.character("eval"), PACKAGE = "RMINC")
    if (length(test) > 1) {
        if (reduce == TRUE) {
            maskV <- mincGetVolume(mask)
            results <- results[maskV > (minmask - 0.5) & maskV < (maxmask + 0.5), ]
        }
        class(results) <- c("mincMultiDim", "matrix")
    }
    else {
        if (reduce == TRUE) {
            maskV <- mincGetVolume(mask)
            results <- results[maskV > (minmask - 0.5) & maskV < (maxmask + 0.5)]
        }
        class(results) <- c("mincSingleDim", "numeric")
    }
    attr(results, "likeVolume") <- filenames[1]
    return(results)
}

