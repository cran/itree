\name{residuals.itree}
\alias{residuals.itree}
\title{
  Residuals From a Fitted itree Object 
}
\usage{
\method{residuals}{itree}(object, type = c("usual", "pearson", "deviance"), ...)
}
\description{
  Method for \code{residuals} for an \code{itree} object.
  This function is used and functions identically to \code{predict.rpart} even for procedures
  unique to itree. Note however that itree is missing the procedures implemented 
  by rpart that are specific to the poisson and exp methods.
}
\arguments{
\item{object}{
  fitted model object of class \code{"itree"}.
}
\item{type}{
  Indicates the type of residual desired.

  For regression or \code{anova} trees all three residual
  definitions reduce to \code{y - fitted}.  This is the residual returned for
  \code{user} method trees as well.

  For classification trees the \code{usual} residuals
  are the misclassification losses L(actual, predicted) where L is the
  loss matrix.  With default losses this residual is
  0/1 for correct/incorrect classification.
  The \code{pearson} residual is 
  (1-fitted)/sqrt(fitted(1-fitted)) and the \code{deviance} residual is 
  sqrt(minus twice logarithm of fitted).
}
\item{\dots}{further arguments passed to or from other methods.}

}
\value{
  Vector of residuals of type \code{type} from a fitted \code{itree} object. 
}
\examples{
#the rpart example:
fit <- itree(skips ~ Opening + Solder + Mask + PadType + Panel,
       data=solder, method='anova')
summary(residuals(fit))
plot(predict(fit),residuals(fit))
}
\keyword{tree}
