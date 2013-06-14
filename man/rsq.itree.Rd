\name{rsq.itree}
\alias{rsq.itree}
\title{
  Plots the Approximate R-Square for the Different Splits 
}
\description{
  Produces 2 plots.  The first plots the r-square (apparent and
  apparent - from cross-validation) versus the number of splits.  
  The second plots the Relative Error(cross-validation) +/- 1-SE from
  cross-validation versus the number of splits. Same as \code{rsq.rpart}
  but does some checking to make sure warnings/error is printed
  when the user attempts to call the function in cases where either
  \code{cptable=NULL} or does not have the correct meaning to make the plots
  useful. Identical to the \code{rpart} function.
}
\usage{
rsq.itree(x)
}
\arguments{
\item{x}{
  fitted model object of class \code{itree}.  This is assumed to be the result
  of some function that produces an object with the same named
  components as that returned by the \code{itree} function. 
}}
\section{Side Effects}{
  Two plots are produced.
}
\note{
  The labels are only appropriate for the \code{"anova"} method. Further
  the cptable from which the r-squared values are taken are not appropriate
  if a penalty has been used.  If this is the case, the method stops.
}
\examples{
#rpart's example:
z.auto <- itree(Mileage ~ Weight, car.test.frame)
rsq.itree(z.auto)
}
\keyword{tree}
