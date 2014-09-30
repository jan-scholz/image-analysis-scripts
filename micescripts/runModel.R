runmodel <- function(table,mask,form,outbase,columns=NA,random='',anova=FALSE, adjust='fdr', verbose=F) {
# formula and random are strings
# write all columns except F-statistic if columns==NA
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

		if (verbose) {
			cat('formula',form,'\n')
			cat('random',random,'\n')
		}

		attach(table)
		cat('running mincLme \n')
		out <- mincLme(filename, form, random, mask=mask)
		detach(table)
	}

	if (adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")) {
		# lme outpout has no df and statType fields
		#q <- mincFDR(out, mask=mask, method="FDR", statType=c('F','t','t','t','t'), df=rep(10,5))
		q <- mincFDR(out, mask=mask, method="FDR")
	}

	colnames <- dimnames(out)[[2]]

	if (is.null(columns)) {
		selectedcolumns <- c()
	} else if (is.na(columns)) {
		start <- 1
		if (colnames[start] == 'F-statistic') start <- 2
		selectedcolumns <- start:length(colnames)
	} else {
		for (x in columns) {
			if(!is.numeric(x) | x>length(colnames)) {
				stop('not a valid column:',x)
			}
		}
		selectedcolumns <- columns
	}

	for (c in selectedcolumns) {
		outname <- paste(outbase,gsub('[()]','',colnames[c]),sep='_')
		mincWriteVolume(out, paste(outname,".mnc",sep=""),     c,clobber=T)
		if (adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")) {
			mincWriteVolume(  q, paste(outname,"_fdr.mnc",sep=""), c,clobber=T)
		}
		cat('saved', outname, '\n')
	}

	con <- file(paste(outbase,".txt",sep=""), open="wt")
	writeLines(paste('# FORMULA: ',toString(form),sep=''), con)
	writeLines(paste('# COLUMNS: ',paste(selectedcolumns,collapse=','),sep=''),con)
	write.csv(attr(q,'thresholds'), con)
	close(con)

	return(out)
}

