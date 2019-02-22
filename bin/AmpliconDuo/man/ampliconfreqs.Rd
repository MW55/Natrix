\name{ampliconfreqs}
\alias{ampliconfreqs}
\docType{data}
\title{
SSU Amplicon Frequencies Of Different Sampling Sites (Short Version)
}
\description{
A data frame with 16 amplicon data sets taken from 8 different sampling sites. A reduced version (frequencies of 2500 amplicons) of the \code{\link{ampliconfreqs.long}} data provided with this package. 
Samples from each sampling site were sequenced twice, corresponding to sets A and B (e.g. column  FU25.A and FU25.B). 
Names of the sampling sites are specified in the data \code{\link{site.f}}.
}
\usage{data(ampliconfreqs)}
\format{
  A data frame with 80903 observations on the following 16 variables.
  \describe{
    \item{\code{FU25.A}}{a numeric vector}
    \item{\code{FU25.B}}{a numeric vector}
    \item{\code{FU28.A}}{a numeric vector}
    \item{\code{FU28.B}}{a numeric vector}
    \item{\code{FU31.A}}{a numeric vector}
    \item{\code{FU31.B}}{a numeric vector}
    \item{\code{FU31.C}}{a numeric vector}
    \item{\code{FU31.D}}{a numeric vector}
    \item{\code{FU34.A}}{a numeric vector}
    \item{\code{FU34.B}}{a numeric vector}
    \item{\code{FU37.A}}{a numeric vector}
    \item{\code{FU37.B}}{a numeric vector}
    \item{\code{UniPond.A}}{a numeric vector}
    \item{\code{UniPond.B}}{a numeric vector}
    \item{\code{BogSoil.A}}{a numeric vector}
    \item{\code{BogSoil.B}}{a numeric vector}
  }
}

\references{
Lange A, Jost S, Heider D, Bock C, Budeus B, Schilling E, Strittmatter A, Boenigk J, Hoffmann D: AmpliconDuo: A Split-Sample Filtering Protocol for High-Throughput Amplicon Sequencing of Microbial Communities. 
\emph{(PLoS One. 2015 Nov 2;10(11))}
}
\examples{
data(ampliconfreqs)
data(site.f)
ampliconduo.a <- ampliconduo(ampliconfreqs[,1:4], sample.names = site.f[1:2])
}
\keyword{datasets}
