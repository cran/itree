\name{estNodeRisk}
\alias{estNodeRisk}
\title{
Estimates risk at a each leaf node in a tree. 
}
\usage{
estNodeRisk(tree.obj,est_observation_loss)
}
\description{
Takes point-wise estimates of loss at each observation used to fit a tree
and aggregates them over the leaf nodes in the tree.
}
\arguments{
	\item{tree.obj}{
	  An itree object whose nodes we wish to have risk estimates for.
	}
	\item{est_observation_loss}{
	A vector of estimated loss for each observation used to fit \code{tree.obj},
	in the same order as \code{tree.obj$y}. This argument is intended to be
	 
	\code{getOOBLoss(tree.obj,mydata,...)$avgOOBloss} 
	but that is up to the user.
	It is up to the USER to ensure that both \code{tree.obj}
	and \code{est_observation_loss} were generated using the same data. Generally 
	they should be generated using the same method and controls as well, but
	this is not required.
	}
}
\value{ List with elements: \cr 
	\code{est.risk} is a vector of length equal to number of leaf nodes in \code{tree.obj}, each element is
	an estimate of out-of-sample risk at a particular node in \code{tree.obj}. This is computed
	by averaging the OOB losses over observations in that node.  The order of \code{est.risk} is
	such that its first element corresponds to the first row that is a leaf node in \code{tree.obj$frame}. 
	est.risk[2] corresponds to the second row that is a leaf node in \code{tree.obj$frame}, and so on. \cr
	
	\code{sd.loss} = sd of the estimated losses for observations belonging to each node. 
	In the same order as \code{est.risk}. 
}
\examples{
require(mlbench); data(BostonHousing)
#fit a tree:
cart <- itree(medv~.,BostonHousing,minsplit=25,minbucket=25,cp=0)

#generate theta-hat values by computing average out-of-bag loss:
\dontrun{
#don't run because of time to do bootstrap sampling...
theta_hats <- getOOBLoss(model_tree.obj=cart.bh,data=bh,nboot=100)

# Then for each leaf we estimate local risk by the mean in-node theta-hat.
lre <- estNodeRisk(tree.obj=cart.bh,est_observation_loss=theta_hats$avgOOBloss)

# to add the lre to the plot:
plot(cart.bh, do_node_re= TRUE, uniform=TRUE)
text(cart.bh, est_node_risk = lre)
}
}
\keyword{tree}