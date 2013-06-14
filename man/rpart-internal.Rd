\name{rpart-internal}
\alias{pred.rpart}
\alias{itree.matrix}
\title{
  Internal Functions 
}
\description{
  Internal functions, only used by packages \pkg{rpart} and \pkg{ipred}.
}
\usage{
pred.rpart(fit, x) 
itree.anova(y, offset, parms, wt)
itree.class(y, offset, parms, wt)
itree.matrix(frame)
itree.class_purity(y, offset, parms, wt)
itree.class_extremes(y, offset, parms, wt)
itree.regression_purity(y, offset, parms, wt)
itree.regression_extremes(y, offset, parms, wt)
}
\arguments{
  \item{fit}{a tree fitted by \code{itree}.}
  \item{x}{a matrix of predictors.}
  \item{y}{the responses.}
  \item{offset}{an offset, or \code{NULL}.}
  \item{parms}{a list of parameters, usually empty.}
  \item{wt}{case weights.}
  \item{frame}{model frame (from call to \code{itree})}
}
\value{
  \code{itree.anova}, \code{itree.class}, \code{itree.exp},
  \code{itree.class_purity},\code{itree.class_extremes},
  \code{itree.regression_purity}, \code{itree.regression_extremes}, 
  \code{itree.poisson} return a list with components
  
  \item{y}{(adjusting for \code{offset} if necessary),}
  \item{parms}{as input,}
  \item{numresp}{the number of responses,}
  \item{summary}{a function to be invoked by \code{\link{summary.itree}},}
  \item{text}{a function to be invoked by \code{\link{text.itree}}.}
}
\keyword{internal}
