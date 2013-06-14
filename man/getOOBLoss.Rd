\name{getOOBLoss}
\alias{getOOBLoss}
\title{
Gets OOB Loss for each observation in a passed data frame.  
}
\usage{
getOOBLoss(model_tree.obj,data,nboot=100,
		sampleFcn = function(idx_vec){sample(idx_vec,replace=TRUE)},
		minsplit, minbucket,lossfcn)
}
\description{
Uses bootstrap sampling to get average loss for each observation. We do this
by averaging out-of-bag loss over a number of runs for each data point. Each in-sample
tree is fit using the same settings as \code{model_tree.obj}.
}
\arguments{
	\item{model_tree.obj}{
	  An \code{itree} object fit using the same variable names as are present in \code{data}.
	  In each bootstrap sample, the parameters of the tree grown to the learning sample
	  are taken to match those of \code{model_tree.obj}.
	}
	\item{data}{
	A data frame -- we construct risk estimates for each row of in this frame by
	looking at average out-of-bag loss.  
	}
	\item{nboot}{ Number of bootrap/cross-validation runs.
	}
	\item{sampleFcn}{ Any function that takes a vector of indices (1,...,N) where
	\code{N=nrow(data)} and returns the indices of a learning sample. The default
	function does bootstrap sampling (size N with replacement). The holdout sample
	is figured out automatically by finding the set difference between 1:N and the
	learning sample.  
	}
	\item{minsplit}{Specifies the \code{minsplit} argument to use for each tree. If blank,
	we get the \code{minsplit} argument from \code{model_tree.obj}. If \code{minsplit}
	is numeric and greater than 1, it represents the number of observations. 
	If \code{minsplit} is numeric and < 1, it is treated as a  
	fraction of N.
	}
	\item{minbucket}{ Same as \code{minsplit} but for the \code{minbucket} argument.
	}
	\item{lossfcn}{ A function that takes two vectors, \code{trueY} and \code{predY}, and outputs
	a vector of losses. If \code{lossfcn} is missing, we default to misclassification rate and
	squared-error-loss for classification and regression respectively. 
	}
}
\value{
	A list with elements: \cr
	
	\code{bagpred} = Average predictions for all N observations over all \code{nboot} runs. For
	classification it is the most common class. \cr 
	
	\code{holdout.predictions} = An \code{N x nboot} matrix with holdout predictions 
	for each run along the columns.  If observation \code{i} is in the learning sample for run j, then
	\code{holdout.predicitions[i,j]} is \code{NA}. \cr
	
	\code{avgOOBloss} = \code{N x 1} vector of average OOB loss. 
}
\examples{
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
}
\keyword{tree}