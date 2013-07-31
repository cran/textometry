\name{textometry-package}
\alias{textometry-package}
\alias{textometry}
\docType{package}

\title{
Statistical exploration of textual corpora
}

\description{
Statistical exploration of textual corpora using several methods
from french 'lexicometrie' and french 'Data Analysis' schools. It includes
methods for exploring irregularity of distribution of lexicon (or other
features) accross parts of texts; multidimensionnal exploration ('Analyse
des correspondances'), etc.
}

\details{
\tabular{ll}{
Package: \tab textometry\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-01-18\cr
License: \tab GPLv3\cr
Depends: \tab R (>= 2.5.0)\cr
Packaged: \tab Thu Mar 12 17:00:01 2009; sloiseau\cr
Built: \tab R 2.7.1; ; 2009-03-12 17:00:03; unix\cr
}

Index:
\preformatted{
specificities    Compute Lexical Specificness of subcorpus 
progression    Draw progression graphic 
}
}

\author{
Sylvain Loiseau, Lise Vaudor, Matthieu Decorde, Lise Vaudor
}

\keyword{ package }

\examples{
data(robespierre);
specificities(robespierre);
}