\name{specificities.probabilities.vector}
\Rdversion{1.1}
\alias{specificities.probabilities.vector}

\title{
  Vector raw hypergeometric probabilities
}
\description{
  Calculate specificity probabilities on vector (call phyper and phyper_right)
}
\usage{
  specificities.probabilities.vector(v_f, v_F, T, t)
}

\arguments{
  \item{v_f}{vector of lexicon ferquencies}
  \item{v_F}{vector of corpus frequencies}
  \item{T}{corpus size}
  \item{t}{sub-corpus size}
}

\value{
  Hypergeometric probabilities. See \code{\link{specificities.lexicon}}.
}

\author{ Matthieu Decorde, Serge Heiden }

