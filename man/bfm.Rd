\name{bfm}
\alias{bfm}
\docType{data}

\title{adverbs frequency from 5 different domains of the BFM database}

\description{
  A lexical table containing frequencies of adverbs from the BFM (Base de Francais
  m\'edi\'eval) database in 5 different domains (literary, historical, didactic, law, religious).
}

\usage{data(bfm)}

\format{
  The format is:
 num [1:2, 1:5] 103000 1370887 23429 413441 15345 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:2] "ADV" "other"
  ..$ : chr [1:5] "literary" "history" "didactic" "juridical" ...
}

\details{
  The last line of the table gives the total frequency of all the other part of speech words in
  each of these domains.
}

\source{
  BFM: \url{https://bfm.ens-lyon.fr}
}

\references{
  BFM - Base de Fran\,cais M\'edi\'eval [En ligne]. Lyon : ENS de Lyon, Laboratoire ICAR, 2012, \url{https://bfm.ens-lyon.fr}.
}

\keyword{datasets}

