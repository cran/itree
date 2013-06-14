\name{post.itree}
\alias{post.itree}
\alias{post}
\title{
  PostScript Presentation Plot of an itree Object 
}
\description{
Generates a PostScript presentation plot of an \code{itree} object. 
Identical to the function in rpart.
}
\usage{
post(tree, \dots)

\method{post}{itree}(tree, title., 
     filename = paste(deparse(substitute(tree)), ".ps", sep = ""), 
     digits = getOption("digits") - 3, pretty = TRUE,
     use.n = TRUE, horizontal = TRUE, \dots)
}
\arguments{
\item{tree}{
fitted model object of class \code{itree}.  This is assumed to be the result 
of some function that produces an object with the same named 
components as that returned by the \code{itree} function. 
}
\item{title.}{
a title which appears at the top of the plot.  By default, the 
name of the \code{itree} endpoint is printed out. 
}
\item{filename}{
ASCII file to contain the output.  By default, the name of the file is 
the name of the object given by \code{itree} (with the suffix \code{.ps} added). If 
\code{filename = ""}, the plot appears on the current graphical device. 
}
\item{digits}{
number of significant digits to include in numerical data. 
}
\item{pretty}{
an integer denoting the extent to which factor levels will be 
abbreviated in the character strings defining the splits; 
(0) signifies no abbreviation of levels.  A \code{NULL} signifies using 
elements of letters to represent the different factor levels.  The 
default (\code{TRUE}) indicates the maximum possible abbreviation. 
}
\item{use.n}{
Logical.  If \code{TRUE} (default), adds to label 
(\#events level1/ \#events level2/etc. for method \code{class}, \code{n} for 
method \code{anova}.
}
\item{horizontal}{
Logical.  If \code{TRUE} (default), plot is horizontal. If \code{FALSE}, plot 
appears as landscape. 
}
\item{\dots}{
other arguments to the \code{postscript} function. 
}
}
\section{Side Effects}{
a plot of \code{itree} is created using the \code{postscript} driver, or
the current device if \code{filename = ""}.
}
\seealso{
\code{\link{plot.itree}},  \code{\link{itree}},  \code{\link{text.itree}},  \code{\link{abbreviate}}   
}
\examples{
#the rpart examples:
z.auto <- itree(Mileage ~ Weight, car.test.frame)
post(z.auto, file = "")   # display tree on active device
   # now construct postscript version on file "pretty.ps"
   # with no title
post(z.auto, file = "pretty.ps", title = " ")
z.hp <- itree(Mileage ~ Weight + HP, car.test.frame)
post(z.hp)
}
\keyword{tree}
