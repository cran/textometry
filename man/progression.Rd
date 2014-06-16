\name{progression}
\alias{progression}

\title{ Draw progression graphic }
\description{
  Draw the progression graphic of matches of CQL queries in a corpus
}
\usage{
  progression(positions, names, colors, styles, widths, corpusname, Xmin, T, 
  	doCumulative, structurepositions, strutnames, graphtitle, bande)
}

\arguments{
  \item{positions}{ Vector containing corpus positions of CQL queries matches. A position is an integer from 0 (begining of corpus) to N (end of corpus) }
  \item{names}{ String vector containing the CQL queries }
  \item{colors}{ Vector containing the line color of each query }
  \item{styles}{ Vector containing the line style of each query }
  \item{widths}{ Vector containing the line width of each query }
  \item{corpusname}{ String: corpus name }
  \item{Xmin}{ Integer: corpus starting position of abscissa values }
  \item{T}{ Integer: size of the corpus }
  \item{doCumulative}{ Boolean: if true draw a cumulative graph, if false draw a density graph }
  \item{structurepositions}{ optional Vector containing the structure positions of the corpus }
  \item{strutnames}{ optional Vector containing the structures labels to display }
  \item{graphtitle}{ String: graph title }
  \item{bande}{ Float: density window size factor }
}

\author{Matthieu Decorde}
