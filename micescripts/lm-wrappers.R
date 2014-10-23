
lmerp <- function(..., nsim=0) {
	library("lme4")
	library("languageR")

	effectname <- names(lmer(...,doFit=F)$fr$fixef)
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

