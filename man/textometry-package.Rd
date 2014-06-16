\name{textometry-package}
\alias{textometry-package}
\alias{textometry}
\docType{package}

\title{
Textual Data Analysis Package used by the TXM Software
}

\description{
Statistical exploration of textual corpora using several methods
from French 'Textometrie' (new name of 'Lexicometrie') and French 'Data Analysis' schools.
It includes methods for exploring irregularity of distribution of lexicon features across
text sets or parts of texts (Specificity analysis); multi-dimensional exploration (Factorial analysis), etc. 
Those methods are used in the TXM software.}

\details{
\tabular{ll}{
Package: \tab textometry\cr
Type: \tab Package\cr
Version: \tab 0.1.3\cr
Date: \tab 2014-06-16\cr
License: \tab GPLv3\cr
Depends: \tab R (>= 1.5.0)\cr
}

Index:
\preformatted{
specificities    Compute Lexical Specificity of subcorpus 
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
