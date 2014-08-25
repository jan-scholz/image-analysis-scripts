# run a command for a certain split mask
# this is run with
#R CMD BATCH --no-restore --no-save --slave '--args tablefile="Data.csv" functionstring="" maskname="foo" maskval=666' ~/bin/micescripts/runformask.R out.Rout && cat out.Rout

# tablefile   CSV file, contains column 'filename' and the ones recquired for functionstring 
# maskname    name of minc mask split up according to integer label
# maskval     the label of the specific submask that will be processed


args=(commandArgs(TRUE))

##args is now a list of character vectors
for(i in 1:length(args)){
	eval(parse(text=args[[i]]))
}


runformask <- function(tablefile, functionstring, maskname, maskval=1) {

	# reduce =TRUE?
	out <- mincApply(t$Absolute_Jacobians,quote( f() ),mask=maskname,maskval=maskval)

#ncol(out)

#outname <- sprintf('%s%04i',basename,maskval)

#	if (column > 0) {
#		selectedcolumns <- column
#	} else {
#		selectedcolumns <- seq(1:length(colnames))
#	}
#
#	for (c in selectedcolumns) {
#		outname <- paste(outbase,gsub('[()]','',colnames[c]),sep='_')
#		mincWriteVolume(out, paste(outname,".mnc",sep=""),     c,clobber=T)
#
#		#q <- mincFDR(out, mask=mask, method="FDR")
#		#mincWriteVolume(q, paste(outname,"_fdr.mnc",sep=""), c,clobber=T)
#
#		cat('saved', outname, '\n')
#	}
#
#	con <- file(paste(outbase,".txt",sep=""), open="wt")
#	writeLines(paste('# FORMULA: ',toString(formula),sep=''), con)
#	#write.csv(attr(q,'thresholds'), con)
#	close(con)



}


runformask(maskname, maskval)

# mincmath -add [<in1.mnc> ...] <out.mnc>

