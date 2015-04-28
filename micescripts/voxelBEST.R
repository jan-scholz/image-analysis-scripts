# library(BEST)
#

voxelBEST <- function(x, data=NULL) {
#voxelBEST <- function(x) {

if (!is.null(data)) d <- data

group         <- 'condition'
ROPEm         <- 0.01
ROPEeff       <- 0.1

numSavedSteps <- 1e3  #1e5
thinSteps     <- 1
burnInSteps   <- 1e3  #1e3

y1 <- x[as.numeric(d[, group])==1]
y2 <- x[as.numeric(d[, group])==2]

# run MCMC
bout <- BESTmcmc(y1, y2, numSavedSteps=numSavedSteps, thinSteps=thinSteps, burnInSteps=burnInSteps, verbose=FALSE)

sout <- summary(bout, ROPEm=c(-ROPEm,ROPEm), ROPEeff=c(-ROPEeff,ROPEeff))

# pre-allocate (faster)
#out <- rep(0,29)
#names(out) <- c("mu1.mean","mu1.mode","mu1.HDIlo","mu1.HDIup","muDiff.mean","muDiff.mode","muDiff.HDIlo","muDiff.HDIup","sigma1.mean","sigma1.mode","sigma1.HDIlo","sigma1.HDIup","sigmaDiff.mean","sigmaDiff.mode","sigmaDiff.HDIlo","sigmaDiff.HDIup","log10nu.mean","log10nu.mode","log10nu.HDIlo","log10nu.HDIup","effSz.mean","effSz.mode","effSz.HDIlo","effSz.HDIup","muDiff.%>compVal","sigmaDiff.%>compVal","effSz.%>compVal","muDiff.%InROPE","effSz.%InROPE")

out <- rep(0,17)
names(out) <- c("mu1.mean","muDiff.mean","muDiff.HDIlo","muDiff.HDIup","sigma1.mean","sigmaDiff.mean","sigmaDiff.HDIlo","sigmaDiff.HDIup","log10nu.mean","effSz.mean","effSz.HDIlo","effSz.HDIup","muDiff.%>compVal","sigmaDiff.%>compVal","effSz.%>compVal","muDiff.%InROPE","effSz.%InROPE")

# assemble output vector
c( out["mu1.mean"]            <- sout["mu1","mean"],
   #out["mu1.mode"]            <- sout["mu1","mode"],
   #out["mu1.HDIlo"]           <- sout["mu1","HDIlo"],
   #out["mu1.HDIup"]           <- sout["mu1","HDIup"],
   out["muDiff.mean"]         <- sout["muDiff","mean"],
   #out["muDiff.mode"]         <- sout["muDiff","mode"],
   out["muDiff.HDIlo"]        <- sout["muDiff","HDIlo"],
   out["muDiff.HDIup"]        <- sout["muDiff","HDIup"],
   out["sigma1.mean"]         <- sout["sigma1","mean"],
   #out["sigma1.mode"]         <- sout["sigma1","mode"],
   #out["sigma1.HDIlo"]        <- sout["sigma1","HDIlo"],
   #out["sigma1.HDIup"]        <- sout["sigma1","HDIup"],
   out["sigmaDiff.mean"]      <- sout["sigmaDiff","mean"],
   #out["sigmaDiff.mode"]      <- sout["sigmaDiff","mode"],
   out["sigmaDiff.HDIlo"]     <- sout["sigmaDiff","HDIlo"],
   out["sigmaDiff.HDIup"]     <- sout["sigmaDiff","HDIup"],
   out["log10nu.mean"]        <- sout["log10nu","mean"],
   #out["log10nu.mode"]        <- sout["log10nu","mode"],
   #out["log10nu.HDIlo"]       <- sout["log10nu","HDIlo"],
   #out["log10nu.HDIup"]       <- sout["log10nu","HDIup"],
   out["effSz.mean"]          <- sout["effSz","mean"],
   #out["effSz.mode"]          <- sout["effSz","mode"],
   out["effSz.HDIlo"]         <- sout["effSz","HDIlo"],
   out["effSz.HDIup"]         <- sout["effSz","HDIup"],
   out["muDiff.%>compVal"]    <- sout["muDiff","%>compVal"],
   out["sigmaDiff.%>compVal"] <- sout["sigmaDiff","%>compVal"],
   out["effSz.%>compVal"]     <- sout["effSz","%>compVal"],
   out["muDiff.%InROPE"]      <- sout["muDiff","%InROPE"],
   out["effSz.%InROPE"]       <- sout["effSz","%InROPE"] )

return(out)
}



#leveneTestForRMINC <- function(x) {
#  fout <- as.numeric(leveneTest(x, gf$Strain, center=mean)[1,2])
#  pout <- as.numeric(leveneTest(x, gf$Strain, center=mean)[1,3])
# return(c(fout, pout))
#}
 

# outLeveneTest <- pMincApply(gf$jacobians, quote(leveneTestForRMINC(x)),  mask="mask.mnc", cores=8)

