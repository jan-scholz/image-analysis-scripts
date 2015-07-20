# model comparison using RMINC's lmer

runModelComp <- function(formList, d, outbase, mask, ncpus=2, bsvoxels=50) {

library(lme4)

out <- vector("list",length(formList))
out_corr <- vector("list",length(formList))
out_p <- vector("list",length(formList))

for (i in 1:length(formList)) {
	cat('########## '); print(formList[[i]])
	out[[i]] <- do.call('mincLmer', list(formula=formList[[i]], data=quote(d), mask=mask, parallel=c("snowfall", ncpus), REML=FALSE))
	out[[i]] <- mincLmerEstimateDF(out[[i]])
	out_corr[[i]] <- mincFDR(out[[i]])

	#cat('i',i,'\n')
	#print(out[[i]])
	#cat('..\n')
	#print(attr(out[[i]],'df'))
	#cat('...\n')

	cnames <- colnames(out[[i]])
	dfs <- attr(out[[i]],'df')

	for (j in 1:ncol(out[[i]])) {
		if (grepl('^tvalue-',cnames[j])) {
			cat('j', j, cnames[j], 'df', dfs[j-length(dfs)], '\n')
			out_p[[i]] <- pt(-abs(out[[i]]), dfs[j-length(dfs)])*2
		}
	}

	for (c in 1:length(cnames)) {
		if (!grepl('Intercept',cnames[c])) {
			mincWriteVolume(out[[i]],      paste(outbase, '_m', i, '_', cnames[c], ".mnc",     sep=''), c, clobber=T, like.filename=mask)
			mincWriteVolume(out_p[[i]],      paste(outbase, '_m', i, '_', cnames[c], "_p.mnc",     sep=''), c, clobber=T, like.filename=mask)
		}
	}

	cnames <- colnames(out_corr[[i]])
	for (c in 1:length(cnames)) {
		if (!grepl('Intercept',cnames[c])) {
			mincWriteVolume(out_corr[[i]], paste(outbase, '_m', i, '_', gsub('qvalue-','',cnames[c]), "_fdr.mnc", sep=''), c, clobber=T, like.filename=mask)
		}
	}
}

#cat(sprintf('%s',Sys.time()),'start mincLogLikRatio\n')
#llr <- do.call('mincLogLikRatio', out)                    do.call with list 'out' very slow for some reason
#cat(sprintf('%s',Sys.time()),'end mincLogLikRatio\n')

llr <- mincLogLikRatio(out[[1]],out[[2]])
llr_corr <- mincFDR(llr, mask=mask)
mincWriteVolume(llr,      paste(outbase, '_llr',     ".mnc", sep=''),    clobber=T)
mincWriteVolume(llr_corr, paste(outbase, '_llr_fdr', ".mnc", sep=''), 1, clobber=T)

# anticonservative bias?
cat('\n\nparametric bootstrap with',bsvoxels,'voxels\n')
llr <- mincLogLikRatioParametricBootstrap(llr, selection = "random", nsims = 500, nvoxels = bsvoxels)
#print(attr(llr, "parametricBootstrap"))
#print(coef(attr(llr$bootstrap,'parametricBootstrapModel')))
#print(confint(attr(llr$bootstrap,'parametricBootstrapModel'),level=c(0.95)))

llr_bscorr <- mincFDR(llr)

out[['bootstrap']] <- llr
out[['bootstrap_corr']] <- llr_bscorr

mincWriteVolume(llr_bscorr, paste(outbase, '_llr_boot_fdr', ".mnc", sep=''), 2, clobber=T)

return(out)

}

# system.time(llr <- mincLogLikRatio2(o[[1]],o[[2]]))
# system.time(llr <- do.call('mincLogLikRatio2', o))

