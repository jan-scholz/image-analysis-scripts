runmodel <- function(table,mask,form,outbase,columns=NA,random='',anova=FALSE, adjust='fdr', colname=T, verbose=F,pvals=F) {
# formula and random are strings
# write all columns except F-statistic if columns==NA
# colname==T : descriptive column names, numbers otherwise
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

		cat('\n#########################################################################\n')
		cat('HINT: make sure variable names do not clash with functions, e.g. "scan()"\n')
		cat('#########################################################################\n\n')

		attach(table)
		cat('running mincLme \n')
		cat('with mask',mask,'\n')
		out <- mincLme(filename, form, random, mask=mask)
		detach(table)
	}

	if (adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")) {
		cat('adjusting\n')
		# lme outpout has no df and statType fields
		#out.q <- mincFDR(out, mask=mask, method="FDR", statType=c('F','t','t','t','t'), df=rep(10,5))
		out.q <- mincFDR(out, mask=mask, method="FDR")
	}

	colnames <- dimnames(out)[[2]]
	if (verbose) {
		cat('\ncolumn names: ', paste(1:length(colnames),colnames, sep=') ', collapse=', '), '\n\n')
	}
	selectedcolumns <- c()

	if (is.null(columns)) {
		if (verbose) cat('no columns selected\n')
		selectedcolumns <- c()
	} else if (is.na(columns)) {
		for (x in 1:length(colnames)) {
			if ( !grepl('(Intercept)',colnames[x]) & !grepl('F-statistic',colnames[x]) & !grepl('R-squared',colnames[x]) & !grepl('beta-',colnames[x]) ) {
				if (verbose) cat(sprintf('selecting column for saving: %2i) %s\n', x, colnames[x]))
				selectedcolumns <- c(selectedcolumns,x)
			}
		}
#		start <- 1
#		if (colnames[start] == 'F-statistic') start <- 2
#		selectedcolumns <- start:length(colnames)
	} else {
		for (x in columns) {
			if(!is.numeric(x) | x>length(colnames)) {
				stop('not a valid column:',x)
			}
		}
		selectedcolumns <- columns
	}
	if (verbose) cat('\n')

	for (c in selectedcolumns) {
		if (colname==T) {
			outname <- paste(outbase,gsub('[()]','',colnames[c]),sep='_')
			#outname <- paste(gsub(' - ','-',outname),sep='_')
			outname <- gsub(' ','_',outname)
		} else {
			outname <- paste(outbase,c,sep='_')
		}

		if (any(is.nan(out))) {
			cat('Warning: NaN detected, setting to 0.\n')
			out[is.nan(out)] <- 0
		}

		if (any(is.infinite(out))) {
			cat('Warning: Inf detected, setting to 0.\n')
			out[is.infinite(out)] <- 0
		}

		if (any(is.na(out))) {
			cat('Warning: NA detected, setting to 0.\n')
			out[is.na(out)] <- 0
		}

		mincWriteVolume(out, paste(outname,".mnc",sep=""), colnames[c], clobber=T)

		if (adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")) {
			if (colnames[c] %in% colnames(out.q)) {
				mincWriteVolume(  out.q, paste(outname,"_fdr.mnc",sep=""), colnames[c], clobber=T, like.filename=mask)
			}
		}

		if (pvals) {
			if (attr(out,'stat-type')[c]=='t') {
				# two-sided p value
				out.p <- pt(-abs(out), attr(out,'df')[[1]][2])*2
				mincWriteVolume(out.p, paste(outname,"_p.mnc",sep=""), colnames[c], clobber=T, like.filename=mask)
			}
		} else {
			cat('incompatible stat type for t-to-p conversion:',colnames[c],'\n')
		}


		cat('saved', outname, '\n\n')
	}

	con <- file(paste(outbase,".txt",sep=""), open="wt")
	writeLines(paste('# FORMULA: ',toString(form),sep=''), con)
	#writeLines(paste('# COLUMNS: ',paste(selectedcolumns,collapse=','),sep=''),con)
	#writeLines(sprintf('# COLUMNS: '), con)
	writeLines(sprintf('# COLUMNS: %s', paste(1:length(colnames),colnames, sep='] ', collapse=',  ')),con)
	write.csv(attr(out.q,'thresholds'), con)
	close(con)

	return(out)
}

