# splits mask into n parts, with ascending integer value
# interleaved starts at 'middle' slice and works outwards (great for intermed. results)

splitmask <- function (mask, basename='', n = 2, interleaved = FALSE) 
{
	if (!file.exists(mask)) stop('Could not find mask: ', mask)
	maskV <- mincGetVolume(mask)
    nVoxels <- sum(maskV > 0.5)

#    if (tinyMask != FALSE) {
#        maskV[maskV > 0.5] <- as.integer(cut(seq_len(nVoxels), tinyMask))
#        maskV[maskV > 1.5] <- 0
#        nVoxels <- sum(maskV > 0.5)
#    }

	if (n < 2) {
		voxellabels <- rep(1,nVoxels)
	} else {
		if (interleaved) {
			interlseq   <- c(seq(n,1,by=-2),seq((n %% 2) + 1,n,by=2))
			voxellabels <- as.numeric(as.character(cut(seq_len(nVoxels), n, labels=interlseq)))
		} else {
			voxellabels <- cut(seq_len(nVoxels), n, labels=FALSE)
		}
	}

	maskV[maskV > 0.5] <- voxellabels

	if (nchar(basename) < 1) {
		maskFilename <- tempfile(pattern="tmpMask", tmpdir=getwd(), fileext=".mnc")
	} else {
		maskFilename <- paste(basename,'.mnc',sep='')
	}

    mincWriteVolume(maskV, maskFilename, clobber=TRUE)
}

# n<-4; c(seq(n,1,by=-2),seq((n %% 2) + 1,n,by=2));

