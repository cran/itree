\name{rpconvert}
\alias{rpconvert}
\title{
  Update an itree object
}
\description{
\code{itree} objects changed in their internal format in order to
accommodate the changes for user-written split functions.  This routine
updates an old object to the new format. Same as \code{rpconvert} 
from the \code{rpart} package.
}
\usage{
rpconvert(x)
}
\arguments{
\item{x}{
an \code{itree} object
}
}
\value{
An updated object
}
\seealso{
\code{itree}
}
\keyword{tree}
