# splits mask into n parts, with ascending integer value

splitmask <- function (mask, basename='', n = 2, order = 'descending') 
{
	if (n < 2) stop('n needs to be at least 2')
	if (!file.exists(mask)) stop('Could not find mask: ', mask)

	library(RMINC)

	maskV <- mincGetVolume(mask)
    nVoxels <- sum(maskV > 0.5)

#    if (tinyMask != FALSE) {
#        maskV[maskV > 0.5] <- as.integer(cut(seq_len(nVoxels), tinyMask))
#        maskV[maskV > 1.5] <- 0
#        nVoxels <- sum(maskV > 0.5)
#    }


	if (order == 'descending') {
		cat('desc\n')
		l <- cut(seq_len(nVoxels), n, labels=FALSE)
		l <- -l+max(l)+1
	} else if (order == 'ascending') {
		cat('asc\n')
		l <- cut(seq_len(nVoxels), n, labels=FALSE)
	} else if (order == 'centrefirst') {
		s <- seq(from=-n+1,to=n-1,by=2)
		s <- abs(s) + as.numeric(sign(s)<0)
		l <- cut(seq_len(nVoxels),n,labels=s)
		l <- as.numeric(as.character(l))
		l[l==0] <- 1
	}

    maskV[maskV > 0.5] <- l

	if (nchar(basename) < 1) {
		maskFilename <- tempfile(pattern="tmpMask", tmpdir=getwd(), fileext=".mnc")
	} else {
		maskFilename <- paste(basename,'.mnc',sep='')
	}


    mincWriteVolume(maskV, maskFilename, clobber=TRUE)
}

