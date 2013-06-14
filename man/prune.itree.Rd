\name{prune.itree}
\alias{prune.itree}
\alias{prune}
\title{
  Cost-complexity Pruning of an itree Object 
}
\description{
Determines a nested sequence of subtrees of the supplied \code{itree} object 
by recursively \code{snipping} off the least important splits, based on the 
complexity parameter (\code{cp}). Identical to \code{prune.rpart}.
}
\usage{
prune(tree, \dots)

\method{prune}{itree}(tree, cp, \dots)
}
\arguments{
\item{tree}{
fitted model object of class \code{itree}.  This is assumed to be the result 
of some function that produces an object with the same named 
components as that returned by the \code{itree} function.
}
\item{cp}{
  Complexity parameter to which the \code{itree} object will be trimmed.
}
\item{\dots}{further arguments passed to or from other methods.}
}
\value{
A new \code{itree} object that is trimmed to the value \code{cp}. 
}
\seealso{
\code{\link{itree}} 
}
\examples{
#the rpart example:
z.auto <- itree(Mileage ~ Weight, car.test.frame)
zp <- prune(z.auto, cp=0.1)
plot(zp) #plot smaller itree object
}
\keyword{tree}
