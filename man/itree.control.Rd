\name{itree.control}
\alias{itree.control}
\title{
  Control for itree Models
}
\description{
Various parameters that control aspects of the \code{itree} fit. 
This is based on \code{rpart.control} but extends it a bit to deal with
the procedures unqiue to \code{itree}.
}
\usage{
itree.control(minsplit = 20, minbucket = round(minsplit/3),
              maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, 
              xval = 10,surrogatestyle = 0, maxdepth = 30, 
              impscale=3,interp_param1=0,interp_param2=0,cp,\dots)
}
\arguments{
  \item{minsplit}{
    the minimum number of observations that must exist in a node in order for 
    a split to be attempted. 
  }
  \item{minbucket}{
    the minimum number of observations in any terminal \code{<leaf>} node.
    If only one of \code{minbucket} or \code{minsplit} is specified,
    the code either sets \code{minsplit} to \code{minbucket*3}
    or \code{minbucket} to \code{minsplit/3}, as appropriate.
  }
  \item{maxcompete}{
    the number of competitor splits retained in the output.  It is useful to
    know not just which split was chosen, but which variable came in second,
    third, etc.
  }
  \item{maxsurrogate}{
    the number of surrogate splits retained in the output.  If this is set to
    zero the compute time will be reduced, since approximately half of the
    computational time (other than setup) is used in the search for surrogate
    splits.
  }
  \item{usesurrogate}{
    how to use surrogates in the splitting process.  \code{0} means
    display only; an observation with a missing value for the primary
    split rule is not sent further down the tree.  \code{1} means use
    surrogates, in order, to split subjects missing the primary variable;
    if all surrogates are missing the observation is not split.  For value
    \code{2} ,if all surrogates are missing, then send the observation in
    the majority direction.  A value of \code{0} corresponds to the action
    of \code{tree}, and \code{2} to the recommendations of Breiman et.al.
  }
  \item{xval}{
    number of cross-validations. This uses the xval method in rpart, which cross-
    validated across a set of \code{cp} values. For one-sided methods \code{cp} is not defined 
    and so passing \code{cp}>0 results in an error.  If a penalty is used, that same
    penalty and penalization constant is used in the cross-validations.
  }
  \item{surrogatestyle}{
    controls the selection of a best surrogate.
    If set to \code{0} (default) the program uses the total number of correct
    classification for a potential surrogate variable,
    if set to \code{1} it uses the percent correct, calculated over the
    non-missing values of the surrogate.
    The first option more severely penalizes covariates with a large number of
    missing values.
  }
  \item{maxdepth}{
    Set the maximum depth of any node of the final tree, with the root
    node counted as depth 0.  Values greater than  30 \code{itree} will
    give nonsense results on 32-bit machines.
  }
  \item{impscale}{
  	How to scale the 'improve' function so that the interpretability penalty can be between 0 and 1.
  	\code{impscale=1} means scale by root node's impurity, \code{=2} means scale by parent node's impurity.
  	\code{=3} means no scaling.  Currently if a penalization method is used and impscale is not 
  	specified, we default to \code{impscale=2}. The default is recommended.
  }
  \item{interp_param1}{
  	First interpretability parameter controlling the tradeoff between
  	reducing loss and introducing less interpretable splits into the tree. Higher values
  	means greater penalty for less interpretable splits.
  }
  \item{interp_param2}{
  	Second interpretability parameter. Currently not used.
  }
    \item{cp}{
    complexity parameter.  Any split that does not decrease the overall
    lack of fit by a factor of \code{cp} is not attempted.  For instance,
    with \code{anova} splitting, this means that the overall Rsquare must
    increase by \code{cp} at each step. The main role of this parameter
    is to save computing time by pruning off splits that are obviously
    not worthwhile. Essentially,the user informs the program that any
    split which does not improve the fit by \code{cp} will likely be
    pruned off by cross-validation, and that hence the program need
    not pursue it. Note that this is not defined for one-sided methods,
    so passing \code{cp=.01} with \code{method="purity"}, for instance, results in an 
    error.  To control the size of purity or extremes trees, use \code{minsplit}
    and/or \code{minbucket}.
  }
  \item{\dots}{
    mop up other arguments.
  }
}
\value{
  A list containing the options.
}
\seealso{
  \code{\link{itree}}
}
\keyword{tree}
