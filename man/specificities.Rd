\name{specificities}
\alias{specificities}

\title{ Calculate Lexical Specificity Score }

\description{
Calculate the specificity - or association or surprise -
score of a word being present \code{f} times or more
in a sub-corpus of \code{t} words given that it appears
a total of \code{F} times in a whole corpus of \code{T} words.
}

\usage{
specificities(lexicaltable, types=NULL, parts=NULL)
}

\arguments{
  \item{lexicaltable}{ a complete lexical table, i.e. a numeric matrix where
  each line represents a word and each column a part of the corpus. Each cell gives the
  frequency of the given word in the corresponding part of the corpus. }

  \item{types}{ list of rows (words) for which the specificity score must be calculated.
  If \code{NULL}, the specificity score is calculated for every row;
  If \code{types} is a character vector, it indicates the row names
  for which the specificity score is to be calculated
  (an error is thrown if \code{lexicaltable} has no row names);
  If \code{types} is an integer vector, it indicates the index of rows
  for which the specificity score is to be calculated.}

  \item{parts}{ list of columns (parts) for which the specificity score must be calculated.
  If \code{NULL}, the specificity index is calculated for every part;
  If \code{parts} is a character vector, it indicates the column names
  for which the specificity score is to be calculated
  (an error is thrown if \code{lexicaltable} has no column names);
  If \code{parts} is an integer vector, it indicates the index of columns
  for which the specificity score is to be calculated.}
}

\value{
Returns a matrix of \code{nrow(lexicaltable) * ncol(lexicaltable)} (the number of
rows and columns may be reduced using \code{types} or \code{parts}), each cell
giving the specificity score.
}

\references{
  Lafon P. (1980) Sur la variabilit\'e de la fr\'e quence des formes dans un corpus,
  Mots, 1, pp. 127--165. \url{https://www.persee.fr/doc/mots_0243-6450_1980_num_1_1_1008}
}

\author{Matthieu Decorde, Serge Heiden, Sylvain Loiseau, Lise Vaudor}

\seealso{ \code{\link{specificities.probabilities}}, \code{\link{specificities.lexicon}} }

\examples{
data(robespierre);
spe <- specificities(robespierre);
string <- paste("The word \%s appears f=\%d times in a sub-corpus of t=\%d words,",
" given a total frequency of F=\%d in the robespierre corpus made",
" of T=\%d words. The corresponding specificity score is \%f", sep="");
print(sprintf(string,
'peuple',
robespierre['peuple','D4'],
colSums(robespierre)['D4'],
rowSums(robespierre)['peuple'],
sum(robespierre),
spe['peuple', 'D4']));
}
