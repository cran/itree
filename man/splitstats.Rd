\name{splitstats}
\alias{splitstats}
\title{
  Stats about a tree's splits.
}
\usage{
splitstats(tree,featlist=NULL)
}
\description{
Computes some quantities about the variables used in the splits of an 
\code{itree} object.  This is new to itree.
}
\arguments{
	\item{tree}{
	  an \code{itree} object.
	}
	\item{featlist}{
	  A list of features to compute statistics of. Typically set to \code{NULL}, in 
	  which case we figure out the features using the \code{itree} object itself.
	}
}
\value{
	A data frame with each row being a feature and each column being a statistic
	about how that feature is used in the tree.
}
\examples{
require(mlbench); data(BostonHousing)
#fit a tree:
cart <- itree(medv~.,BostonHousing,minsplit=25,minbucket=25,cp=0)
splitstats(cart)
}
\keyword{tree}