anatModel <- function(data,formula,atlas,anova=F,combined=F,adjust='fdr') {
	# anatModel() retrieves volumes for each atlas label and returns the p-value
	#             for the model specified by the formula
	#
	# copyright 2013 jan.scholz@phenogenomics.ca @ MICe, Toronto, Canada
	# version 2013-03-13
	#
	# data          data frame with regressors; contains filename column
	# formula       formula
	# atlas         filename of MINC file with labels
	# anova         if true, use ANOVA
	# combined      combine left and right
	# adjust        correction method, e.g. 'none', 'fdr', see p.adjust for more
	#
	# value         returns array of p-values

	method = "jacobians"

	library(RMINC)
	library(reshape)

	if ( ! 'filename' %in% colnames(data) ) {
		stop('data does not contain column \'filename\'')
	}

	# get volume data
	vols <- anatGetAll(data$filename, atlas=atlas, method=method)
	if (combined) {
		vols <- anatCombineStructures(vols, method=method)
	}

	data <- melt(cbind(data, vols), id=seq(ncol(data)))
	labelnames <- levels(data$variable)

	# correct formula if wrong left side
	formula <- formula(gsub('^.*~ *','value ~ ', deparse(formula)))

	# initalize output array
	if (anova) {
		effectnames <- rownames(anova(lm(formula, data=subset(data, variable==labelnames[1]))))
	} else {
		effectnames <- names(coef(lm(formula, data=subset(data, variable==labelnames[1]))))
	}
	p <- array(NA, dim=c(ncol(vols), length(effectnames)), dimnames=list(labelnames,effectnames))

	# populate array with p values
	for (s in rownames(p)) {
		if (anova) {
			p[s,] <- anova(lm(formula, data=subset(data,variable==s)))$`Pr(>F)`
		} else {
			#c[s] <- coef(lm(formula, data=subset(data,variable==s)))
			p[s,] <- coef(summary(lm(formula, data=subset(data,variable==s))))[,'Pr(>|t|)']
		}
	}

	# adjusted p values
	if (adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")) {
		for (c in colnames(p)) {
			p[,c] <- p.adjust(p[,c], adjust)
		}
	} else {
		adjust <- 'none'
	}

	attr(p, "p.adjust.method") <- adjust
	attr(p, "formula") <- deparse(formula)
	attr(p, "anova") <- anova

	return(p)
}

