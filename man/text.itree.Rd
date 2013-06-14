\name{text.itree}
\alias{text.itree}
\title{
  Place Text on a Dendrogram 
}
\description{
  Labels the current plot of the tree dendrogram with text.
  Extends \code{text.rpart} to also print \code{itree}'s node risk estimate
  along with the fitted values.
}
\usage{
\method{text}{itree}(x, splits=TRUE, label, FUN=text, all=FALSE,
     pretty=NULL, digits=getOption("digits") - 3, use.n=FALSE,
     fancy=FALSE, fwidth=0.8, fheight=0.8, est_node_risk = NULL,use_sd=FALSE, \dots)
}
\arguments{
\item{x}{
  fitted model object of class \code{itree}.  This is assumed to be the result
  of some function that produces an object with the same named
  components as that returned by the \code{itree} function.
}
\item{splits}{
  logical flag.  If \code{TRUE} (default), then the splits in the tree are
  labeled with the criterion for the split.
}
\item{label}{
% a column name of \code{x$frame};  values of this
% will label the nodes.  For the \code{"class"} method,
% \code{label="yval"} results in the factor levels
% being used, \code{"yprob"} results in the probability of the winning factor
% level being used, and \sQuote{specific yval level} results in the probability
% of that factor level.
For compatibility with \code{itree2}, ignored in this version (with
a warning).
}
\item{FUN}{
  the name of a labeling function, e.g. \code{text}.  
}
\item{all}{
  Logical.  If \code{TRUE}, all nodes are labeled, otherwise just
  terminal nodes.
}
\item{pretty}{
  an integer denoting the extent to which factor levels in split labels
  will be abbreviated.  A value of (0) signifies no abbreviation.  A
  \code{NULL}, the default, signifies using elements of letters to represent
  the different factor levels.
}
\item{digits}{
  number of significant digits to include in numerical labels.
}
\item{use.n}{
  Logical.  If \code{TRUE}, adds to label 
  (\#events level1/ \#events level2/etc. for \code{class},
  \code{n} for \code{anova}, and \#events/n for \code{poisson} and \code{exp}).
}
\item{fancy}{
  Logical.  If \code{TRUE}, nodes are represented by ellipses (interior nodes)
  and rectangles (leaves) and labeled by \code{yval}.  The edges connecting the
  nodes are labeled by left and right splits.
}
\item{fwidth}{
  Relates to option \code{fancy} and the width of the ellipses and rectangles.
  If \code{fwidth < 1} then it is a scaling factor (default = 0.8).  If 
  \code{fwidth > 1} then it represents the number of character widths 
  (for current graphical device) to use.
}
\item{fheight}{
  Relates to option \code{fancy} and the height of the ellipses and rectangles.
  If \code{fheight <1} then it is a scaling factor (default = 0.8).  If
  \code{fheight > 1} then it represents the number of character heights 
  (for current graphical device) to use.
}
\item{est_node_risk}{
If not \code{NULL}, an object of class \code{estNodeRisk} giving the risk estimates
to print at each node.
}
\item{use_sd}{
If \code{TRUE}, prints sd of the average losses for each node. Default is \code{FALSE}.
}
\item{\dots}{
  Graphical parameters may also be supplied as arguments to this
  function (see \code{par}).  As labels often extend outside the plot
  region it can be helpful to specify \code{xpd = TRUE}.
}}
\section{Side Effects}{
  the current plot of a tree dendrogram is labeled.
}
\seealso{
  \code{\link{text}}, \code{\link{plot.itree}}, \code{\link{itree}},
  \code{\link{post.itree}}, \code{\link{abbreviate}} 
}
\examples{
# from rpart:
freen.tr <- itree(y ~ ., freeny)
plot(freen.tr)
text(freen.tr, use.n=TRUE, all=TRUE)

###unique to itree, plotting node risk:
require(mlbench); data(BostonHousing)
#fit a tree:
cart <- itree(medv~.,BostonHousing,minsplit=25,minbucket=25,cp=0)

#generate theta-hat values by computing average out-of-bag loss:
\dontrun{
#don't run it because of time to do the bootstrap...
theta_hats <- getOOBLoss(model_tree.obj=cart.bh,data=bh,nboot=100)

# Then for each leaf we estimate local risk by the mean in-node theta-hat.
lre <- estNodeRisk(tree.obj=cart.bh,est_observation_loss=theta_hats$avgOOBloss)

# to add the lre to the plot:
plot(cart.bh, do_node_re= TRUE, uniform=TRUE)
text(cart.bh, est_node_risk = lre)
}
}
\keyword{tree}
