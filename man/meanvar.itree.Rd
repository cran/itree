\name{meanvar.itree}
\alias{meanvar}
\alias{meanvar.itree}
\title{
  Mean-Variance Plot for an itree Object 
}
\description{
Creates a plot on the current graphics device of the deviance of the 
node divided by the number of observations at the node.  Also returns 
the node number. Same as \code{meanvar.rpart}, but implements checks so that
a warning/error is printed when called on an object where deviance
doesn't have the appropriate meaning.
}
\usage{
meanvar(tree, \dots)

\method{meanvar}{itree}(tree, xlab="ave(y)", ylab="ave(deviance)", \dots)
}
\arguments{
\item{tree}{
fitted model object of class \code{itree}.  This is assumed to be the result 
of some function that produces an object with the same named 
components as that returned by the \code{itree} function. 
}
\item{xlab}{
x-axis label for the plot. 
}
\item{ylab}{
y-axis label for the plot. 
}
\item{\dots}{
additional graphical parameters may be supplied as arguments to this function. 
}
}
\value{
an invisible list containing the following vectors is returned. 

\item{x}{
fitted value at terminal nodes (\code{yval}). 
}
\item{y}{
deviance of node divided by number of observations at node. 
}
\item{label}{
node number. 
}
}
\section{Side Effects}{
a plot is put on the current graphics device. 
}
\seealso{
 \code{\link{plot.itree}}.   
}
\examples{
#the rpart package's example:
z.auto <- itree(Mileage ~ Weight, car.test.frame)
meanvar(z.auto, log='xy')
}
\keyword{tree}
