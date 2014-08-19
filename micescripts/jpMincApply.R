# TODO: catch ctrl-c and clean up temporary masks

jpMincApply <- function (filenames, function.string, mask = NULL, cores = 4, 
    tinyMask = FALSE, method = "snowfall") 
{
    if (is.null(mask)) {
        stop("err, need a mask. Sorry. Will fix it soon")
    }
    else {
        maskV <- mincGetVolume(mask)
    }
    nVoxels <- sum(maskV > 0.5)
    if (tinyMask != FALSE) {
        maskV[maskV > 0.5] <- as.integer(cut(seq_len(nVoxels), tinyMask))
        maskV[maskV > 1.5] <- 0
        nVoxels <- sum(maskV > 0.5)
    }
    maskV[maskV > 0.5] <- as.integer(cut(seq_len(nVoxels), cores))

	maskFilename <- tempfile(fileext='.mnc')
    #maskFilename <- paste("pmincApplyTmpMask-", Sys.getpid(), ".mnc", sep = "")
    mincWriteVolume(maskV, maskFilename, clobber = TRUE)
    pout <- list()
    test <- eval(function.string)
    if (method == "local") {
        stop("Lovely code ... that generates inconsistent results because something somewhere is not thread safe ...")
        library(multicore)
        library(doMC)
        library(foreach)
        registerDoMC(cores)
        pout <- foreach(i = 1:cores) %dopar% {
            mincApply(filenames, function.string, mask = maskFilename, 
                maskval = i)
        }
    }
    else if (method == "sge") {
        stop("implementation of sge method completely broken ...")
        i <- 4
        pout <- list()
        pids <- sge.submit(mincApply, filenames, function.string, 
            mask = maskFilename, maskval = i, packages = c("RMINC"))
        status <- sge.job.status(pids$pid)
        while (status != 0) {
            Sys.sleep(4)
            status <- sge.job.status(pids$pid)
        }
        pout[[i]] <- sge.list.get.result(pids)
    }
    else if (method == "snowfall") {
        wrapper <- function(i) {
            return(mincApply(filenames, function.string, mask = maskFilename, 
                maskval = i, reduce = TRUE))
        }
        if (is.null(cores)) {
            cores <- length(sfSocketHosts())
        }
        pout <- sfLapply(1:cores, wrapper)
    }
    else {
        stop("unknown execution method")
    }
    if (length(test) > 1) {
        output <- matrix(0, nrow = length(maskV), ncol = length(test))
        class(output) <- class(pout[[1]])
        attr(output, "likeVolume") <- attr(pout[[1]], "likeVolume")
    }
    else {
        output <- maskV
    }
    for (i in 1:cores) {
        if (length(test) > 1) {
            output[maskV == i, ] <- pout[[i]]
        }
        else {
            output[maskV == i] <- pout[[i]]
        }
    }
    unlink(maskFilename)
    return(output)
}

