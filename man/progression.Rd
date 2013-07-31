\name{progression}
\alias{progression}

\title{ Draw progression graphic }
\description{
  Draw the progression of matches of queries in the corpus
}
\usage{
  progression(positions, names, colors, styles, widths, corpusname, Xmin, T, 
  	doCumulative, structurepositions, strutnames, graphtitle, bande)
}

\arguments{
  \item{positions}{ A vector containing the matching positions of queries }
  \item{names}{ A string vector containing the queries }
  \item{colors}{ A vector containing the colors per query }
  \item{styles}{ A vector containing the styles per query }
  \item{widths}{ A vector containing the widths per query }
  \item{corpusname}{ The corpus name as String}
  \item{Xmin}{ Where starting the abscissa values}
  \item{T}{ Integer: The size of the corpus }
  \item{doCumulative}{ Boolean: if true draw do cumulative graph, if false draw a density graph }
  \item{structurepositions}{ optional Vector containing the structures position in the corpus }
  \item{strutnames}{ optional Vector containing the structures label to display}
  \item{graphtitle}{ String: graph title}
  \item{bande}{ Float: multiply the density window.}
}

\value{
  a progression plot
}

\author{Matthieu Decorde}