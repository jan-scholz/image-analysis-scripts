runLmer <- function(table,mask,formula,outbase,column=0,cores=4) {

	doLmer <<- function(x) {
		l <- lmer(as.formula(formula), data=table)
		# return t values for every column
		return(coef(summary(l))[,3])
	}

	# initalize prallel environment
	library(snowfall)
	sfInit(parallel=TRUE, cpus=cores)
	sfLibrary(lme4)
	sfClusterEval(setwd(getwd()))
	sfExport('table','doLmer','mask')

	out <- pMincApply(table$filename,quote(doLmer(x)),mask=mask,cores=cores)
	#out <- pMincApply(table$filename,quote(doLmer(x)),mask=mask,cores=cores,tinyMask=600)
	sfStop()


	# get column names by running on a single voxel
	ts <- mincGetVoxel(table$filename, 0, 0, 0)
	tmp <- doLmer(ts)
	colnames <- names(tmp)

	if (column > 0) {
		selectedcolumns <- column
	} else {
		selectedcolumns <- seq(1:length(colnames))
	}

	for (c in selectedcolumns) {
		outname <- paste(outbase,gsub('[()]','',colnames[c]),sep='_')
		mincWriteVolume(out, paste(outname,".mnc",sep=""),     c,clobber=T)

		#q <- mincFDR(out, mask=mask, method="FDR")
		#mincWriteVolume(q, paste(outname,"_fdr.mnc",sep=""), c,clobber=T)

		cat('saved', outname, '\n')
	}

	con <- file(paste(outbase,".txt",sep=""), open="wt")
	writeLines(paste('# FORMULA: ',toString(formula),sep=''), con)
	#write.csv(attr(q,'thresholds'), con)
	close(con)

	return(out)
}

