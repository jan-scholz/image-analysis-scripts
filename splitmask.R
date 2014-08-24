# splits mask into n parts, with ascending integer value

splitmask <- function (mask, basename='', n = 2) 
{
	if (!file.exists(mask)) stop('Could not find mask: ', mask)
	maskV <- mincGetVolume(mask)
    nVoxels <- sum(maskV > 0.5)

#    if (tinyMask != FALSE) {
#        maskV[maskV > 0.5] <- as.integer(cut(seq_len(nVoxels), tinyMask))
#        maskV[maskV > 1.5] <- 0
#        nVoxels <- sum(maskV > 0.5)
#    }
    maskV[maskV > 0.5] <- cut(seq_len(nVoxels), n, labels=FALSE)

	if (nchar(basename) < 1) {
		maskFilename <- tempfile(pattern="tmpMask", tmpdir=getwd(), fileext=".mnc")
	} else {
		maskFilename <- paste(basename,'.mnc',sep='')
	}


    mincWriteVolume(maskV, maskFilename, clobber=TRUE)
}

