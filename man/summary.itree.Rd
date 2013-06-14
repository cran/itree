\name{summary.itree}
\alias{summary.itree}
\title{
  Summarize a Fitted itree Object 
}
\description{
  Returns a detailed listing of a fitted \code{itree} object.
  Same code and returned value as \code{rpart} but with additional 
  info regarding penalties and any new parameters specific to \code{itree} objects.
}
\usage{
\method{summary}{itree}(object, cp=0, digits=getOption("digits"), file, \dots)
}
\arguments{
\item{object}{
  fitted model object of class \code{itree}.  This is assumed to be the result
  of some function that produces an object with the same named
  components as that returned by the \code{itree} function.
}
\item{digits}{
    Number of significant digits to be used in the result.
}
\item{cp}{
  trim nodes with a complexity of less than \code{cp} from the listing.
}
\item{file}{
  write the output to a given file name. (Full listings of a tree are
  often quite long).
}
\item{\dots}{
  arguments to be passed to or from other methods.
}}
\details{
  This function is a method for the generic function summary for class
  \code{"itree"}.  It can be invoked by calling \code{summary}
  for an object of the appropriate class, or directly by calling
  \code{summary.itree} regardless of the class of the object.
}
\seealso{
  \code{\link{summary}}, \code{\link{itree.object}}, \code{\link{printcp}}.
}
\examples{
#rpart's example:
z.auto <- itree(Mileage ~ Weight, car.test.frame)
summary(z.auto)
}
\keyword{tree}
