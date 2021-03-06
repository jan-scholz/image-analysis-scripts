# run a command for a certain split mask
# functionstring does not seem to work, at least not always ...
# 
# run permuted stats with lmert
# R --no-save --no-restore -e "source('~/bin/micescripts/runformask2.R'); runformask(outbase='dti_faCmeanFA_lmert/out0003', tablefile='table_all_local_fas.csv', functionstring='lmert(x~sex+group+meanFA+(1|side)+(1|id),data=datatable)', mask='splitmask80.mnc', maskval=1, libs=c('lme4'), verbose=TRUE)"

runformask <- function(outbase, tablefile, functionstring='', mask, maskval=1, columns=NA, libs=c(), colname=T, permute=list(idvar='',filenamecols=c()), verbose=FALSE) {
	# colname==T : descriptive column names, numbers otherwise

	NotFiniteReplacement <- 0
	source('~/bin/micescripts/lm-wrappers.R')

	if (!file.exists(tablefile)) stop('ERROR: Could not find table file: ',tablefile)
	if (!file.exists(mask)) stop('ERROR: Could not find mask file: ',mask)

	library(RMINC)
	for (l in libs) if (l != '') library(l ,character.only = TRUE)

	v <- mincGetVolume(mask)
	#if (sum(v>(maskval-0.5) & v<(maskval+0.5))<1) stop('ERROR: mask does not contain mask value: ',maskval)
	v <- array(v,dim=rev(minc.dimensions.sizes(mask)))   # reverse order here (works when mincinfo yields dim order: z,y,xspace)
	firstmaskvoxel <- which((v>(maskval-0.5) & v<(maskval+0.5)), arr.ind=T)[1,]-1
	rm(v)

	datatable <<- read.csv(tablefile)

	# sample
	if (nchar(permute$idvar) & length(permute$filenamecols)) {
		if (verbose) cat('##### permuting input data\n')
		datatable <- permuteData(datatable, idvar=permute$idvar, filenamecols=permute$filenamecols)   # permute filename-related data
		cat('##### permuted EV\n')
		print(datatable[,c('oldid','filename')])
		#datatable <- sampleLmeData(datatable, idvar='id')
	}

	x <<- mincGetVoxel(datatable$filename, firstmaskvoxel)
	singlevoxelmodel <- eval(parse(text=functionstring))
	if (verbose) cat('##### finished running single voxel model\n')

	colnames <- names(singlevoxelmodel)
	selectedcolumns <- c()

	if (is.null(columns)) {
		if (verbose) cat('##### no output\n')
	} else if (is.na(columns)) {
		if (verbose) cat('##### excluding "(Intercept)" from output\n')
		for (x in 1:length(colnames)) {
			if ( !grepl('(Intercept)',colnames[x]) ) {
				selectedcolumns <- c(selectedcolumns,x)
			}
		}
	} else {
		for (x in columns) {
			if(!is.numeric(x) | x>length(colnames)) {
				stop('not a valid column:',x)
			}
		}
		selectedcolumns <- columns
	}

	if (verbose) {
		for (x in selectedcolumns) {
				if (verbose) cat(' #### adding column to output:', colnames[x], '\n')
		}
	}

	# MAIN ###################################################################
	f <- paste('quote(',functionstring,')',sep='')
	starttime <- Sys.time()
	out <- mincApply(datatable$filename,eval(parse(text=f)),mask=mask,maskval=maskval)
	endtime <- Sys.time()
    ##########################################################################

	for (c in selectedcolumns) {

		if (colname==T) {
			outname <- paste(outbase,gsub('[()]','',colnames[c]),sep='_')
			#outname <- paste(gsub(' - ','-',outname),sep='_')
			outname <- gsub(' ','_',outname)
		} else {
			outname <- paste(outbase,c,sep='_')
		}

		if (sum(!is.finite(out)) > 0) {
			warning('replacing non-numbers with', NotFiniteReplacement)
			out <- replace(out,!is.finite(out),NotFiniteReplacement)
		}
		mincWriteVolume(out, paste(outname,".mnc",sep=""),     c,clobber=T)
		cat('saved', outname, '\n')
	}

	con <- file(paste(outbase,".txt",sep=""), open="wt")
	writeLines(paste('# FUNCTIONSTRING: ',toString(functionstring),sep=''), con)
	writeLines(paste('# COLUMNS: ',paste(selectedcolumns,collapse=','),sep=''),con)
	writeLines(sprintf('# START: %s\n# END:   %s\n# TIME (MIN): %0.2f',starttime, endtime, difftime(endtime,starttime,units='mins')),con)
	close(con)
}

###############################################################################
# for TFCE
sampleLmeData <- function(dataframe, idvar = "PERSONID") {
	uniqueIDs  <- unique(get(idvar, dataframe))
	sampledIDs <- data.frame(id=sample(uniqueIDs, replace=T))
	newdata    <- merge(dataframe, sampledIDs, by.x=idvar, by.y="id")
	return(newdata)
}


