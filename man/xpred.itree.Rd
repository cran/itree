\name{xpred.itree}
\alias{xpred.itree}
\title{
  Return Cross-Validated Predictions 
}
\description{
  Gives the predicted values for an \code{itree} fit, under
  cross validation, for a set of complexity parameter values.
  Similar to \code{xpred.rpart} but will stop and print an error message
  when given a fit for which \code{cp} is not defined.  If \code{fit} has a penalty,
  then all the cross-validations use the same penalty and penalization 
  constant found in the \code{fit} object.
}
\usage{
xpred.itree(fit, xval=10, cp)
}
\arguments{
\item{fit}{
a \code{itree} object. 
}
\item{xval}{
  number of cross-validation groups.  This may also be an explicit list of
  integers that define the cross-validation groups.
}
\item{cp}{
  the desired list of complexity values.  By default it is taken from the
  \code{cptable} component of the fit.
}}
\value{
  A matrix with one row for each observation and one column for each complexity
  value.
}
\details{
  From rpart:\cr
  Complexity penalties are actually ranges, not values.  If the
  \code{cp} values found in the table were \eqn{.36}, \eqn{.28},
  and \eqn{.13}, for instance, this means that the first row of the
  table holds for all complexity penalties in the range \eqn{[.36, 1]},
  the second row for \code{cp} in the range \eqn{[.28, .36)} and
  the third row for \eqn{[.13,.28)}.  By default, the geometric mean
  of each interval is used for cross validation.
}
\seealso{
\code{\link{itree}} 
}
\examples{
#rpart's example:
fit <- itree(Mileage ~ Weight, car.test.frame)
xmat <- xpred.itree(fit)
xerr <- (xmat - car.test.frame$Mileage)^2
apply(xerr, 2, sum)   # cross-validated error estimate

# approx same result as rel. error from printcp(fit)
apply(xerr, 2, sum)/var(car.test.frame$Mileage) 
printcp(fit)
}
\keyword{tree}
