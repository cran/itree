\name{itree.object}
\alias{itree.object}
\title{
  Recursive Partitioning and Regression Trees Object 
}
\description{
These are objects representing fitted \code{itree} trees. 
}
\section{Structure}{
The following components must be included in a legitimate \code{itree} object. 
Of these, only the \code{where} component has the same length as  
the data used to fit the \code{itree} object.  The requirements here are the same
as those in rpart except itree objects have a \code{penalty} parameter.
}
\value{

\item{frame}{
data frame with one row for each node in the tree.
The \code{row.names} of \code{frame} contain the (unique) node numbers that
follow a binary ordering indexed by node depth.
Elements of \code{frame} include \code{var}, a factor giving
the variable used in the split at each node
(leaf nodes are denoted by the string \code{<leaf>}), \code{n}, the size of each node,
\code{wt}, the sum of case weights for the node,
\code{dev}, the deviance of each node,
\code{yval}, the fitted value of the response at each node,
and \code{splits}, a two column matrix of left and right split labels
for each node. All of these are the same as for an \code{itree} object. \cr

For classification problems, information about total counts (or weights,
if weights are unequal) appear in the \code{wt.classX} column where
the integer \code{X} ranges from 1 to the number of classes.  Similarly, 
the \code{wt.frac.classX} is the weight of class X in the node divided by
the total weight in the node. \code{nodewt} is the total weight of 
all observations in this node as fraction of the entire dataset. This naming
convention is different from rpart's.\cr

Also included in the frame are \code{complexity}, the
complexity parameter at which this split will collapse, \code{ncompete},
the number of competitor splits retained, and \code{nsurrogate}, the
number of surrogate splits retained. Note that \code{complexity} values are
dependent on any penalty method and penalization constant used. 
}
\item{where}{
integer vector, the same length as the number of observations in the root node, 
containing the row number of \code{frame} corresponding to the leaf node 
that each observation falls into. 
}
\item{splits}{
a numeric matrix describing the splits.  The row label is the name of the split
variable, and columns are \code{count}, the number of observations sent left
or right by the split (for competitor splits this is the number that
would have been sent left or right had this split been used, for surrogate
splits it is the number missing the primary split variable which were decided
using this surrogate), \code{ncat}, the number of categories or levels for the
variable (\code{+/-1} for a continuous variable), \code{improve}, which is the improvement
in deviance given by this split, or, for surrogates, the concordance of the
surrogate with the primary, and \code{split}, the numeric split point.  
The last column \code{adj} gives the adjusted concordance for surrogate splits.
For
a factor, the \code{split} column contains the row number of the csplit matrix.
For a continuous variable, the sign of \code{ncat} determines whether the
subset \code{x < cutpoint} or \code{x > cutpoint} is sent to the left.
}
\item{csplit}{
this will be present only if one of the split variables is a factor. There
is one row for each such split, and column \code{i = 1} if this level of the
factor goes to the left, \code{3} if it goes to the right, and 2 if that level
is not present at this node of the tree.
For an ordered categorical variable all levels are marked as \code{R/L}, 
including levels that are not present.
}
\item{method}{
the method used to grow the tree. 
}
\item{penalty}{
the penalty function for splitting on a specific variable at a specific node given
the variables used in the branch leading to this node.
}
\item{cptable}{
the table of optimal prunings based on a complexity parameter. \code{NULL} for 
extremes and purity methods.
}
\item{terms}{
an object of mode \code{expression} and class \code{term} summarizing the formula. 
Used by various methods, but typically not of direct relevance to users. 
}
\item{call}{
an image of the call that produced the object, but with the arguments 
all named and with the actual formula included as the formula argument. 
To re-evaluate the call, say \code{update(tree)}. 

Optional components include the matrix of predictors (\code{x}) and the 
response variable (\code{y}) used to construct the \code{itree} object. 
}
}
\seealso{
\code{\link{itree}}.
}
\keyword{tree}
\keyword{methods}
