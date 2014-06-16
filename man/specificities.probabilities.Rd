\name{specificities.probabilities}
\Rdversion{1.1}
\alias{specificities.probabilities}

\title{
Calculate specificity probabilities
}

\description{
Utility function computing specificity probabilities for the \code{\link{specificities}} function.
}

\usage{
specificities.probabilities(lexicaltable, types = NULL, parts = NULL)
}

\value{
Returns a matrix of signed specificity probabilities (between -1.0 and 1.0).
By convention:
\item{sign}{The sign indicates if the observed frequency is lower (minus) or higher (plus) than the mode of the specificity model}
\item{.Machine$double.xmin limit}{-10.0 and 10.0 values are used to hold the sign when the zero/.Machine$double.xmin boundary line has been crossed (the \code{\link{phyper}} function always returns 0.0)}
}

\arguments{
  \item{lexicaltable}{ see \code{\link{specificities}} }

  \item{types}{ see \code{\link{specificities}} }

  \item{parts}{ see \code{\link{specificities}} }
}

\seealso{
 see \code{\link{specificities}}.
}
