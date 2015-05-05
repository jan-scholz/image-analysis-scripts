# run a command for a certain split mask
# functionstring does not seem to work, at least not always ...

runformask <- function(outbase, tablefile, functionstring='', mask, maskval=1, columns=NA, libs=c(), functions=c(), colname=T, verbose=FALSE) {
	# colname==T : descriptive column names, numbers otherwise

	NotFiniteReplacement <- 0
	source('/micehome/jscholz/bin/micescripts/lm-wrappers.R')

	if (!file.exists(tablefile)) stop('ERROR: Could not find table file: ',tablefile)
	if (!file.exists(mask)) stop('ERROR: Could not find mask file: ',mask)

	library(RMINC)
	for (l in libs) if (l != '') library(l ,character.only = TRUE)
	for (f in functions) if (f != '') source(f)

	v <- mincGetVolume(mask)
	#if (sum(v>(maskval-0.5) & v<(maskval+0.5))<1) stop('ERROR: mask does not contain mask value: ',maskval)
	v <- array(v,dim=rev(minc.dimensions.sizes(mask)))   # reverse order here (works when mincinfo yields dim order: z,y,xspace)
	firstmaskvoxel <- which((v>(maskval-0.5) & v<(maskval+0.5)), arr.ind=T)[1,]-1
	rm(v)

	datatable <<- read.csv(tablefile)
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


