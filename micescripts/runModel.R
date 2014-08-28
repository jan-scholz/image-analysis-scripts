runmodel <- function(table,mask,form,outbase,columns=NULL,random='',anova=FALSE) {
# formula and random are strings
	library(RMINC)

	if (!'filename' %in% colnames(table)) warning('missing column: filename')

	if (random=='' & anova==FALSE) {
		cat('running mincLm \n')
		out <- mincLm(as.formula(form), data=table, mask=mask)
	} else if (random=='' & anova==TRUE) {
		cat('running mincAnova \n')
		out <- mincAnova(as.formula(form), data=table, mask=mask)
	} else if (!random=='' & anova==FALSE) {
		library(nlme)
		form <- gsub('^ *filename *~ *','x ~ ',form)
		if (grepl('[()]',random)) warning('no brackets allowed in random effect')

		cat('form',form,'\n')
		cat('rand',rand,'\n')

		attach(table)
		cat('running mincLme \n')
		out <- mincLme(filename, form, random, mask=mask)
	}

	q <- mincFDR(out, mask=mask, method="FDR")
	colnames <- dimnames(out)[[2]]

	if (!is.null(columns)) {
		selectedcolumns <- columns
	} else {
		start <- 1
		if (colnames[start] == 'F-statistic') start <- 2
		selectedcolumns <- start:length(colnames)
	}

	for (c in selectedcolumns) {
		outname <- paste(outbase,gsub('[()]','',colnames[c]),sep='_')
		mincWriteVolume(out, paste(outname,".mnc",sep=""),     c,clobber=T)
		mincWriteVolume(  q, paste(outname,"_fdr.mnc",sep=""), c,clobber=T)
		cat('saved', outname, '\n')
	}

	con <- file(paste(outbase,".txt",sep=""), open="wt")
	writeLines(paste('# FORMULA: ',toString(form),sep=''), con)
	writeLines(paste('# COLUMNS: ',paste(selectedcolumns,collapse=','),sep=''),con)
	write.csv(attr(q,'thresholds'), con)
	close(con)

	return(out)
}

