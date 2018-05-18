\name{site.f}
\alias{site.f}
\docType{data}
\title{
Sampling Site Short Names
}
\description{
Character vector with the names of the sampling sites, corresponding to
the names used to denote amplicon frequencies in 
the \code{\link{ampliconfreqs}} data.
}
\usage{data(site.f)}
\format{
  The format is:
 chr [1:8] "FU25" "FU28" "FU31.1" "FU31.2" "FU34" "FU37" "UniPond" "BogSoil"
}

\references{
Boenigk J, Heider D, Jost S, Lange A, Budeus B, Schilling E,
Strittmatter A, Hoffmann D: A high-throughput amplicon sequencing and
analysis protocol for comparative analyses of microbial communities
\emph{(submitted)}
}

\examples{
data(site.f)
data(ampliconfreqs)
ampliconduo(ampliconfreqs[,1:6], sample.names = site.f[1:3])

}
\keyword{datasets}
