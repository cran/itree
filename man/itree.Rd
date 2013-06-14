\name{itree}
\alias{itree}
\title{
  Recursive Partitioning and Regression Trees a la rpart,
  but with some extensions targeted at growing interpretable/parsimonious
  trees.
}
\description{
  Fit a \code{itree} model.
}
\usage{
itree(formula, data, weights, subset, na.action = na.itree,
	 method, penalty= "none", model = FALSE, 
	 x = FALSE, y = TRUE, parms, control, cost, \dots)
      
}
\arguments{
  \item{formula}{a \link{formula}, with a response but no interaction terms.}
  
  \item{data}{an optional data frame in which to interpret the variables
    named in the formula.} 
  
  \item{weights}{optional case weights.}
  
  \item{subset}{optional expression saying that only a subset of the
    rows of the data should be used in the fit.}
  
  \item{na.action}{the default action deletes all observations for which
    \code{y} is missing, but keeps those in which one or more predictors
    are missing.}
  
  \item{method}{one of \code{"anova"}, \code{"class"},
    \code{"extremes"}, \code{"purity"}, \code{"class_extremes"}, 
    \code{"class_purity"}, \code{"regression_extremes"}, or 
    \code{"regression_purity"}. The purity and extremes methods
    are new to itree. Unlike rpart, itree does not currently support  
    \code{method="poisson"} or \code{method="exp"}. If \code{method} is missing 
    then the routine tries to make an intelligent guess -- the default is the
    CART methodology, as in rpart. If \code{y} is a factor then 
    \code{method="class"} is assumed, otherwise \code{method="anova"} is assumed.
    Passing a factor with \code{method="purity"} is equivalent to passing
    \code{method="class_purity"}, and similarly for extremes/regression. 
    It is wisest to specify the method directly, especially as more
    criteria may added to the function in future.\cr \cr
    As in rpart, \code{method} can be a list of functions named \code{init}, \code{split} and \code{eval}.
    See the rpart documentation for how this works.
    }
    
  \item{penalty}{one of \code{"none"}, \code{"newvar"} or \code{"ema"}. The penalty
  for splitting a particular node on a specified predictor given the predictors already used in the
   branch leading to this node. Default is "none" which corresponds to CART. \code{"newvar"}
  penalizes predictors not used in the branch leading to the current node. \code{"ema"} implements an
  exponential moving average style penalty whereby recently used variables are favored.
  }
  
  \item{model}{if logical: keep a copy of the model frame in the result?
    If the input value for \code{model} is a model frame (likely from an
    earlier call to the \code{itree} function), then this frame is used
    rather than constructing new data.}
    
  \item{x}{keep a copy of the \code{x} matrix in the result.}
  
  \item{y}{keep a copy of the dependent variable in the result.  If
    missing and \code{model} is supplied this defaults to \code{FALSE}.}
  
  \item{parms}{optional parameters for the splitting function.\cr
    Anova splitting has no parameters.\cr
    For classification splitting, the list can contain any of:
    the vector of prior probabilities (component \code{prior}), the loss matrix
    (component \code{loss}) or the splitting index (component
    \code{split}).  The priors must be positive and sum to 1.  The loss
    matrix must have zeros on the diagonal and positive off-diagonal
    elements.  The splitting index can be \code{gini} or
    \code{information}.  The default priors are proportional to the data
    counts, the losses default to 1, and the split defaults to
    \code{gini}. \cr
    For the regression extremes method, \code{parms=1} or 
    \code{parms=-1 }specifies whether
    we are looking for high or low means respectively (see Buja & Lee). Default is \code{1} 
    for high means.\cr
    For classification extremes, parms is a list specificying the class of interest -- see the
    examples for syntax.
    }
  
  \item{control}{a list of options that control details of the
    \code{itree} algorithm, similar to \code{rpart.control}.  
    See \code{\link{itree.control}}.}
  
  \item{cost}{a vector of non-negative costs, one for each variable in
    the model. Defaults to one for all variables.  These are scalings to
    be applied when considering splits, so the improvement on splitting
    on a variable is divided by its cost in deciding which split to
    choose. Note that costs are not currently supported by the extremes or purity
    methods.}
  
  \item{\dots}{arguments to \code{\link{itree.control}} may also be
    specified in the call to \code{itree}.  They are checked against the
    list of valid arguments.}
}

\details{
  itree is based on the code of rpart, but with some extensions targeted at growing interpretable/parsimonious
  trees.  Bug reports and the like should be directed to this package's maintainer -- not rpart's.
}

\value{
  An object of class \code{itree}.  See \code{\link{itree.object}}.
}

\references{
  Breiman, Friedman, Olshen, and Stone. (1984)
  \emph{Classification and Regression Trees.}
  
  Buja, Andreas and Lee, Yung-Seop (2001). \emph{Data Mining Criteria for Tree-Based Regression and
Classification}, Proceedings of KDD 2001, 27-36.

  Wadsworth.
}

\seealso{
  \code{\link{itree.control}}, \code{\link{itree.object}},
  \code{\link{summary.itree}}, \code{\link{print.itree}}
}

\examples{
#CART (same as rpart):
fit <- itree(Kyphosis ~ Age + Number + Start, data=kyphosis)
fit2 <- itree(Kyphosis ~ Age + Number + Start, data=kyphosis,
              parms=list(prior=c(.65,.35), split='information'))
fit3 <- itree(Kyphosis ~ Age + Number + Start, data=kyphosis,
              control=itree.control(cp=.05))
par(mfrow=c(1,2), xpd=NA) # otherwise on some devices the text is clipped
plot(fit)
text(fit, use.n=TRUE)
plot(fit2)
text(fit2, use.n=TRUE)

#### new to itree:
#same example, but using one-sided extremes:
fit.ext <- itree(Kyphosis ~ Age + Number + Start, data=kyphosis,method="extremes",
				parms=list(classOfInterest="absent"))
#we see buckets with every y="absent":
plot(fit.ext); text(fit.ext,use.n=TRUE) 


library(mlbench); data(BostonHousing)

#one sided purity:
fit4 <- itree(medv~.,BostonHousing,method="purity",minbucket=25)

#low means tree:
fit5 <- itree(medv~.,BostonHousing,method="extremes",parms=-1,minbucket=25)

#new variable penalty:
fit6 <- itree(medv~.,BostonHousing,penalty="newvar",interp_param1=.2)

#ema penalty
fit7 <- itree(medv~.,BostonHousing,penalty="ema",interp_param1=.1)

#one-sided-purity + new variable penalty:
fit8 <- itree(medv~.,BostonHousing,method="purity",penalty="newvar",interp_param1=.2)

#one-sided extremes for classification must specify a "class of interest"
data(PimaIndiansDiabetes)
levels(PimaIndiansDiabetes$diabetes)  
fit9.a <- itree(diabetes~.,PimaIndiansDiabetes,minbucket=50,
                 method="extremes",parms=list(classOfInterest="neg"))
                 
plot(fit9.a); text(fit9.a)

#can also pass the index of the class of interest in levels().
fit9.b <- itree(diabetes~.,PimaIndiansDiabetes,minbucket=50,
                 method="extremes",parms=list(classOfInterest=1))
# so fit9.a = fit9.b
}


\keyword{tree}
