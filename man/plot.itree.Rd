\name{plot.itree}
\alias{plot.itree}
\title{
Plot an itree Object
}
\description{
Plots an \code{itree} object on the current graphics device. This is
based on the plotting function from rpart but modifies/extends
it in various ways to deal with some of \code{itree}'s capabilities.
}
\usage{
\method{plot}{itree}(x, uniform=FALSE, branch=1, compress=FALSE, nspace,
     margin=0, minbranch=.3, highlight.color="black",do_node_re=FALSE, \dots)
}
\arguments{
\item{x}{
  a fitted object of class \code{itree}, containing a classification,
  regression, or rate tree.
}
\item{uniform}{
if \code{TRUE}, uniform vertical spacing of the nodes is used; this may be
less cluttered when fitting a large plot onto a page.
The default is to use a non-uniform spacing proportional to the 
error in the fit.
}
\item{branch}{
controls the shape of the branches from parent to child node.
Any number from 0 to 1 is allowed.  A value of 1 gives square
shouldered branches, a value of 0 give V shaped branches, 
with other values being intermediate.
}
\item{compress}{
if \code{FALSE}, the leaf nodes will be at the horizontal plot coordinates of
\code{1:nleaves}.  
If \code{TRUE}, the routine attempts a more compact arrangement of
the tree.
The compaction algorithm assumes \code{uniform=TRUE}; surprisingly, the result
is usually an improvement even when that is not the case.
}
\item{nspace}{
the amount of extra space between a node with children and
a leaf, as compared to the minimal space between leaves.  
Applies to compressed trees only.  The default is the value of
\code{branch}.
}
\item{margin}{
an extra fraction of white space to leave around the borders of the tree.
(Long labels sometimes get cut off by the default computation).
}
\item{minbranch}{
set the minimum length for a branch to \code{minbranch} times the average
branch length.  This parameter is ignored if \code{uniform=TRUE}.
Sometimes a split will give very little improvement, or even (in
the classification case) no improvement at all.  
A tree with branch lengths strictly proportional to improvement
leaves no room to squeeze in node labels.
}
\item{highlight.color}{
	If the \code{itree} object to be plotted is one-sided, setting highlight.color="RED",
	for example, will highlight the path to the node for which the split was chosen.
	For instance if method="purity", the highlights outline the pure nodes.
}
\item{do_node_re}{
	Set to \code{TRUE} if you wish to call \code{text.itree()} with risk estimates for 
	each leaf node (default is \code{FALSE}). 	This matters in the plot function 
	because space must be left to accommodate writing
	the risk estimates for each terminal node.
}
\item{\dots}{
  arguments to be passed to or from other methods.
}}
\value{
The coordinates of the nodes are returned as a list, with
components \code{x} and \code{y}.
}
\section{Side Effects}{
An unlabeled plot is produced on the current graphics device.
}
\details{
This function is a method for the generic function \code{plot}, for objects
of class \code{itree}. The y-coordinate of the top node of the tree will always be 1.
}
\seealso{
\code{\link{itree}}, \code{\link{text.itree}}
}
\examples{
#the rpart example:
fit <- itree(Price ~ Mileage + Type + Country, cu.summary)
plot(fit, compress=TRUE)
text(fit, use.n=TRUE)

### new to itree, plotting node risk:
require(mlbench); data(BostonHousing)
#fit a tree:
cart <- itree(medv~.,BostonHousing,minsplit=25,minbucket=25,cp=0)

#generate theta-hat values by computing average out-of-bag loss:
\dontrun{
theta_hats <- getOOBLoss(model_tree.obj=cart.bh,data=bh,nboot=100)

# Then for each leaf we estimate local risk by the mean in-node theta-hat.
lre <- estNodeRisk(tree.obj=cart.bh,est_observation_loss=theta_hats$avgOOBloss)

# to add the lre to the plot:
plot(cart.bh, do_node_re= TRUE, uniform=TRUE)
text(cart.bh, est_node_risk = lre)
}

#plot using highlighting for one-sided methods:
purity.tree <- itree(medv~.,BostonHousing,minsplit=25,minbucket=25,cp=0,method="purity")
plot(purity.tree,highlight.color="blue")
text(purity.tree)

}
\keyword{tree}
