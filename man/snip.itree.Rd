\name{snip.itree}
\alias{snip.itree}
\title{
  Snip Subtrees of an \code{itree} Object 
}
\description{
  Creates a "snipped" \code{itree} object, containing the nodes that remain
  after selected subtrees have been snipped off.  The user can snip
  nodes using the toss argument, or interactively by clicking the mouse
  button on specified nodes within the graphics window. Same as \code{snip.rpart}.
}
\usage{
snip.itree(x, toss)
}
\arguments{
\item{x}{
  fitted model object of class \code{itree}.  This is assumed to be the result
  of some function that produces an object with the same named
  components as that returned by the \code{itree} function.
}
\item{toss}{
  an integer vector containing indices (node numbers) of all subtrees to
  be snipped off.  If missing, user selects branches to snip off as
  described below.
}}
\value{
  A \code{itree} object containing the nodes that remain after specified or
  selected subtrees have been snipped off.  
}
\details{
  A dendrogram of \code{itree} is expected to be visible on the graphics
  device, and a graphics input device (e.g., a mouse) is required.  Clicking
  (the selection button) on a node displays the node number, sample
  size, response y-value, and Error (dev). Clicking a second time on the
  same node snips that subtree off and visually erases the subtree.
  This process may be repeated an number of times.  Warnings result from
  selecting the root or leaf nodes.  Clicking the exit button will stop
  the snipping process and return the resulting \code{itree} object.

  See the documentation for the specific graphics device for details on
  graphical input techniques.
}
\section{Warning}{
  Visually erasing the plot is done by over-plotting with the background
  colour.  This will do nothing if the background is transparent
  (often true for screen devices).
}
\seealso{
  \code{\link{plot.itree}}
}
\examples{
## dataset not in R
## same as rpart:
\dontrun{
z.survey <- itree(market.survey) #grow the itree object
plot(z.survey) #plot the tree
z.survey2 <- snip.itree(z.survey,toss=2) #trim subtree at node 2
plot(z.survey2) #plot new tree

# can also interactively select the node using the mouse in the
# graphics window

}
#works with itree fits as well:
z.auto <- itree(Mileage ~ Weight, method="purity",car.test.frame)
#snip off the long branch:
z.auto2 <- snip.itree(z.auto,toss=3)
}
\keyword{tree}
