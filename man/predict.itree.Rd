\name{predict.itree}
\alias{predict.itree}
\title{
  Predictions from a Fitted itree Object 
}
\description{
  Returns a vector of predicted responses from a fitted \code{itree} object. 
  This function is used and functions identically to \code{predict.rpart} even for procedures
  unique to \code{itree}. However \code{itree} does not currently support the \code{type="matrix"}
  option.
}
\usage{
\method{predict}{itree}(object, newdata = list(),
       type = c("vector", "prob", "class", "matrix"),
       na.action = na.pass, \dots)
}
\arguments{
  \item{object}{
    fitted model object of class \code{itree}.  This is assumed to be
    the result of some function that produces an object with the same named 
    components as that returned by the \code{itree} function. 
  }
  \item{newdata}{
    data frame containing the values at which predictions are required. 
    The predictors referred to in the right side of
    \code{formula(object)} must be present by name in \code{newdata}.
    If missing, the fitted values are returned. 
  }
  \item{type}{
    character string denoting the type of predicted value returned.  If
    the \code{itree} object is a classification tree, then the default is to
    return \code{prob} predictions, a matrix whose columns are the
    probability of the first, second, etc. class.  (This agrees with the
    default behavior of \code{\link[tree]{tree}}).  Otherwise, a vector
    result is returned.
  }
  \item{na.action}{a function to determine what should be done with
    missing values in \code{newdata}.  The default is to pass them down
    the tree using surrogates in the way selected when the model was built.
    Other possibilities are \code{na.omit} and \code{na.fail}.
  }
  \item{\dots}{
    further arguments passed to or from other methods.
  }
}
\value{
  A new object is obtained by
  dropping \code{newdata} down the object.  For factor predictors, if an
  observation contains a level not used to grow the tree, it is left at
  the deepest possible node and \code{frame$yval} at the node is the
  prediction.

  If \code{type="vector"}:\cr
  vector of predicted responses.
  For regression trees this is the mean response at the node, for Poisson
  trees it is the estimated response rate, and for classification trees
  it is the predicted class (as a number).
  
  If \code{type="prob"}:\cr
  (for a classification tree) a matrix of class probabilities.

  If \code{type="matrix"}:\cr
  For regression trees, this is the mean response, and for classification 
  trees it is the concatenation of the predicted class, the class counts at that node in the fitted tree, and
  the class probabilities.  

  If \code{type="class"}:\cr
  (for a classification tree) a factor of classifications based on the
  responses.
}
\details{
  This function is a method for the generic function predict for class
  \code{itree}. It can be invoked by calling \code{predict} for an object
  of the appropriate class, or directly by calling \code{predict.itree}
  regardless of the class of the object.
}
\seealso{
\code{\link{predict}},  \code{\link{itree.object}}   
}
\examples{
#verbatim from rpart examples:
z.auto <- itree(Mileage ~ Weight, car.test.frame)
predict(z.auto)

fit <- itree(Kyphosis ~ Age + Number + Start, data=kyphosis)
predict(fit, type="prob")   # class probabilities (default)
predict(fit, type="vector") # level numbers
predict(fit, type="class")  # factor
predict(fit, type="matrix") # level number, class frequencies, probabilities

sub <- c(sample(1:50, 25), sample(51:100, 25), sample(101:150, 25))
fit <- itree(Species ~ ., data=iris, subset=sub)
fit
table(predict(fit, iris[-sub,], type="class"), iris[-sub, "Species"])
}
\keyword{tree}
