# wrapper functions for lmer, lme, etc.
# yield appropiately labeled output array for all derived measure, e.g. estimate, t, p


###############################################################################
# for TFCE
permuteData <- function(dataframe, idvar="id", filenamecols = c("filename")) {

	# create image related data frame, all columns incl. 'filename' related to image data, e.g. mean
	filenamedata <- dataframe[,colnames(dataframe) %in% c(filenamecols,idvar)]

	# non-image related data frame (i.e. remaining columns) containing only id-unique columns
	remainingdata <- unique(dataframe[,!colnames(dataframe) %in% filenamecols])
	remainingdata$oldid <- remainingdata[,idvar]

	# shuffle id's
	remainingdata[,idvar] <- factor(sample(levels(remainingdata[,idvar]),replace=F))

	newdata <- merge(remainingdata,filenamedata,by=idvar)

	return(newdata)
}


###############################################################################
# lmer fast version, no p-value, just t-values, needs runformask's libs=c('lme4')
lmert <- function(formula, data) {
	#library("lme4")
	tryCatch({
		l <- lmer(formula, data)
		out <- attributes(summary(l))$coef[,'t value']
	},
	error=function(e) {
		warning(sprintf('lmerp: %s',e))
		out <- NULL
	},
	finally=return(out)
	) # tryCatch
}


###############################################################################
# lmer with optional p value calculation
lmerp <- function(..., nsim=0) {
	library("lme4")
	library("languageR")

	effectname <- names(fixef(lmer(...)))
	#effectname <- names(lmer(...,doFit=F)$fr$fixef)
#effectname <- names(lmer(x~group+(1|side)+(1|id),data=t,doFit=F)$fr$fixef)
	neffects <- length(effectname)

	statname       <- c('Estimate','t value','Pr(>|t|)', 'pMCMC')
	statname_short <- c(    '-est',     '-t',      '-p','-pMCMC')
	default        <- c(        0,        0,         1,       1)

	nstats <- 2	
	if (nsim > 0) nstats <- 4

	# initialize output with default values and right dimensions
	out <- NULL
	for (i in 1:neffects) {
		for (j in 1:nstats) {
			out[paste(effectname[i],statname_short[j],sep='')] <- default[j]
		}
# put error default 0 in here !!!!!!!!!!!!!
	}

	tryCatch({
		l <- summary(lmer(...))
#l <- summary(lmer(x~group+(1|side)+(1|id),data=t))
		coefs <- attr(l,'coefs')

		# get p value with MCMC
		if (nsim > 0) {
			coefsmcmc <- pvals.fnc(l,addPlot=F,ndigits=6,nsim=nsim)          # 5e5
			coefs <- cbind(coefs,sapply(coefsmcmc$fixed,as.numeric)[,2:6])   # first column 'Estimate' removed
		}

		# save to output
		for (i in 1:neffects) {
			err <- 1
			for (j in 1:nstats) {
				if (is.finite(coefs[effectname[i],statname[j]])) {
					out[paste(effectname[i],statname_short[j],sep='')] <- coefs[effectname[i],statname[j]]
					err <- err*1
				} else {
					err <- err*0
				}
			}
			out[paste(effectname[i],'-err',sep='')] <- (err+0.5)
		}
	},
	error=function(e) { warning(sprintf('lmerp: %s',e)); },
	finally=return(out)
	) # tryCatch
}


###############################################################################
lmep <- function(...) {
	ErrorReplacement <- c(0,0,0,1)
	library("nlme")
	tryCatch({
		l <- summary(lme(...))
		out <- array()
		for (r in rownames(l$tTable)) {
			out[paste(r,'-est',sep='')] <- l$tTable[r,'Value']
			out[paste(r,'-df',sep='')]  <- l$tTable[r,'DF']
			out[paste(r,'-t',sep='')]   <- l$tTable[r,'t-value']
			out[paste(r,'-p',sep='')]   <- l$tTable[r,'p-value']
		}
		# out[-1] first element is NA for whatever reason
		out <- out[-1]
		out <- as.numeric(!is.finite(out))*ErrorReplacement
		return(out)
	},
	error=function(e) { return(ErrorReplacement); }
	) # tryCatch
}


###############################################################################
lmp <- function(...) {
	ErrorReplacement <- c(0,0,1)
	tryCatch({
		l <- summary(lm(...))
		out <- array()
		for (r in rownames(l$coef)) {
			out[paste(r,'-est',sep='')] <- l$coef[r,'Estimate']
			out[paste(r,'-t',sep='')]   <- l$coef[r,'t value']
			out[paste(r,'-p',sep='')]   <- l$coef[r,'Pr(>|t|)']
		}
		# out[-1] first element is NA for whatever reason
		out <- out[-1]
		out <- as.numeric(!is.finite(out))*ErrorReplacement
		return(out)
	},
	error=function(e) { return(ErrorReplacement); }
	) # tryCatch
}




################################################################################
## run lmert multiple times on sampled data frame (with replacement), keep t-values
#lmert_permut <- function(formula, data, nsim=1, idvar = "id") {
#	library("lme4")
#
#	effectname <- names(lmer(formula,data,doFit=F)$fr$fixef)
#	out <- NULL
#
#	for (i in 1:nsim) {
#		data_sample <- sampleLmeData(data, idvar)
#		out <- c(out,lmert(formula, data_sample))  # rbind discards error NULL
#	}
#
#	return(out)
#}


###############################################################################
## for TFCE
#sampleLmeDatasplit <- function(dataframe, splitvar, idvar = "PERSONID") {
#	uniqueIDs  <- unique(get(idvar, dataframe))
#	sampledIDs <- data.frame(id=sample(uniqueIDs, replace=T))
#
#	splitdata <- list()
#	newdata <- NULL
#	for (v in levels(dataframe[,splitvar])) {
#		splitdata[[v]] <- dataframe[dataframe[,splitvar]==v,]
#		newdata <- rbind(newdata,merge(dataframe, sampledIDs, by.x=idvar, by.y="id"))
#	}
#	return(newdata)
#}




###############################################################################

#cd /projects/mice/jscholz/rot/reg_dti+flip/stats
# R --no-save --no-restore -e "source('~/bin/micescripts/runformask.R'); runformask(outbase='dti_faCmeanFA_lmerp/out0019', tablefile='table_all_local_fas.csv', functionstring='lmerp(x~sex+group+meanFA+(1|side)+(1|id),data=datatable,nsim=1e1)', mask='splitmask50.mnc', maskval=19, verbose=TRUE)"
# library(RMINC)
# library(lme4)
# t <- read.csv('table_all_local_fas.csv')
# x <- mincGetVoxel(t$filename,c(50,58,26))
# formula <- x~sex+group+meanFA+(1|side)+(1|id)
# data <- t
# library(ggplot2)
# ggplot(t,aes(group,x,colour=side)) + geom_boxplot()
# l <- lmer(x~group+(1|side)+(1|id),data=t)


###############################################################################
#permuted <- data.frame(id=levels(dataframe[,idvar]),group=sample(levels(dataframe[,permvar]),size=nlevels(dataframe[,idvar]),replace=T))
#newdata    <- merge(dataframe[!colnames(dataframe) %in% permvar], permuted, by.x=idvar, by.y="id")

#sampleLmeData <- function(dataframe, idvar = "PERSONID") {
#	uniqueIDs  <- unique(get(idvar, dataframe))
#	sampledIDs <- data.frame(id=sample(uniqueIDs, replace=T))
#	newdata    <- merge(dataframe, sampledIDs, by.x=idvar, by.y="id")
#	return(newdata)
#}
#


###############################################################################
# for TFCE
#permuteData <- function(dataframe, idvar = "id", permvar = "filename") {
#	# create table with one permvar entry per idvar level (by building contingency table idvar vs permvar)
#	tmp <- data.frame(table(dataframe[,idvar],dataframe[,permvar]))
#	tmp <- tmp[tmp$Freq>0,1:2]  # select idvar, pervar columns
#	colnames(tmp) <- c(idvar,permvar)
#
#	# permute permvar
#	tmp[,permvar] <- sample(tmp[,permvar],replace=F)
#
#	# optionally sample tmp$idvar here
#	#N <- nrow(tmp)
#	#tmp[sample(tmp$id,N,replace=T),]
#
#	newdata <- merge(dataframe[!colnames(dataframe) %in% permvar],tmp, by=idvar)
#	return(newdata)
#}
#
