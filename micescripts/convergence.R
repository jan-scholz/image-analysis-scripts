# calculate convergence time of two crossing time courses
# copyright 2013-02-12, jan.scholz@phenogenomics.ca @ MICe
#
# b <- boot(data=d, statistic=convergence.time, R=1e2, strata=d$group, delta=0.01, epsilon=0.001, parallel='multicore', ncpus=4)

convergence.time <- function(data,indices,delta,epsilon,df) {
	# data: contains time, group, value
	# delta: temporal sampling
	# epsilon: group time courses have converged if abs(difference) < epsilon
	# df: degrees of freedom of spline fit
	# returns: c(convergence_time, convergence_boolean)

	d <- data[indices,]
	l <- lm(value ~ group * ns(time, df=df), data=d)

	tminmax <- range(d$time)
	ngroups <- nlevels(d$group)

	# predict on more fine-grained time interval
	  time <- rep(seq(tminmax[1],tminmax[2],by=delta),ngroups)
	 group <- rep(levels(d$group),each=length(time)/ngroups)
	     p <- data.frame(time,group)
	p$pred <- predict(l,p)

	# difference between and convergence time
	p$diff <- ddply(p,.(time),summarize,diff=diff(pred))$diff
	subthresh <- p[abs(p$diff)<epsilon,'time']

	if (length(subthresh) == 0) {
		return(c(runif(1,tminmax[1],tminmax[2]), 0))
	} else {
		return(c(min(subthresh),                 1))
	}
}


run.convergence.time <- function(x,data,delta=0.5,epsilon=0.01,df=2,R=10,parallel='no',ncpus=1) {
	#data contains: group, time 

	reqcols <- c('id','time','group')
	if (length(intersect(reqcols, names(data))) < 3) stop('data needs to contain columns: id, time, group')

	data$value <- x

	b <- boot(data=data, statistic=convergence.time, R=R, strata=data$group, delta=delta, epsilon=epsilon, df=df, parallel=parallel, ncpus=ncpus)

	bhit <- b$t[b$t[,2]==1,1]
	if (length(bhit) < 1) bhit <- c(-20,20)

	# do not need to return b$t0[1], because all the same values
	return( c( tc=mean(b$t[,1]),
tc_hit=mean(bhit), 
diff(quantile(b$t[,1],c(.1,.9))),
diff(quantile(bhit,c(.1,.9))),
hits=mean(b$t[,2])) )

}


