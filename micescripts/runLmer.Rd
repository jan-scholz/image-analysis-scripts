\name{runLmer}
\alias{runLmer}
\title{wrapper of mincApply on lmer}
\description{
  Runs \code{lmer} on every voxel using mincApply.
}
\usage{
runLmer(table,mask,formula,outbase,column=NA,cores=4)
}
\arguments{
  \item{table}{a data.frame containing a column with paths to Jacobian determinants}
  \item{mask}{mask file name}
  \item{formula}{lmer formula as a string}
  \item{outbase}{basename of output files}
  \item{cores}{number of cores to use}
}
\seealso{
  \code{\link{mincApply}}.
}
\examples{
## one example
prunmodel(table,mask,'filename ~ group * time (1|id)','results')
}
\keyword{lmer}

