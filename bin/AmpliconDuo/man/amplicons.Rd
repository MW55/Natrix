\name{amplicons}
\alias{amplicons}
\docType{data}
\title{
Example List of Ampliconduo Data Frames
}
\description{
A list of ampliconduo data frames that was generated calling the
\code{\link{ampliconduo}} function with its default parameters and the 
\code{\link{ampliconfreqs.long}} data provided with this package as input. 
}
\usage{data(amplicons)}

\references{
Boenigk J, Heider D, Jost S, Lange A, Budeus B, Schilling E,
Strittmatter A, Hoffmann D: A high-throughput amplicon sequencing and
analysis protocol for comparative analyses of microbial communities
\emph{(submitted)}
}

\examples{
data(amplicons)
plotAmpliconduo.set(amplicons, nrow = 3)
}
\keyword{datasets}
