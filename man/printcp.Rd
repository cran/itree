\name{printcp}
\alias{printcp}
\title{
  Displays CP table for Fitted itree Object 
}
\description{
  Displays the \code{cp} table for fitted \code{itree} object.
  Note that cp is not defined for \code{method="purity"} or \code{"extremes"}.
  Otherwise identical to rpart's \code{printcp} function. 
}
\usage{
printcp(x, digits=getOption("digits") - 2)
}
\arguments{
\item{x}{
fitted model object of class \code{itree}.  This is assumed to be the result 
of some function that produces an object with the same named 
components as that returned by the \code{itree} function. 
}
\item{digits}{
  the number of digits of numbers to print.
}}
\details{
Prints a table of optimal prunings based on a complexity parameter. Note that
cp is not defined for one-sided methods. 
}
\seealso{
\code{\link{summary.itree}},  \code{\link{itree.object}}   
}
\examples{
#rpart's examples:

z.auto <- itree(Mileage ~ Weight, car.test.frame)
printcp(z.auto)
\dontrun{
Regression tree:
itree(formula = Mileage ~ Weight, data = car.test.frame)

Variables actually used in tree construction:
[1] Weight

Root node error: 1354.6/60 = 22.576

        CP nsplit rel error  xerror     xstd 
1 0.595349      0   1.00000 1.03436 0.178526
2 0.134528      1   0.40465 0.60508 0.105217
3 0.012828      2   0.27012 0.45153 0.083330
4 0.010000      3   0.25729 0.44826 0.076998
}}
\keyword{tree}
