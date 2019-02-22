\name{AmpliconDuo-package}
\alias{AmpliconDuo-package}
\alias{AmpliconDuo}
\docType{package}
\title{
Statistical Analysis Of Amplicon Data Of The Same Sample To Identify Artefacts
}
\description{

Increasingly powerful techniques for high-throughput sequencing open the
possibility to comprehensively characterize microbial communities,
including rare species. However, a still unresolved issue are the
substantial error rates in the experimental process generating these
sequences. To overcome these limitations we propose an approach, where
each sample is split and the same amplification and sequencing protocol
is applied to both halves. This procedure should allow to detect likely PCR and sequencing artifacts, 
and true rare species by comparison of the results of both parts.

The \pkg{AmpliconDuo} package, whereas ampliconduo from here on refers to the two amplicon data sets of a split sample,
is intended to help interpret the obtained amplicon frequency distribution across split samples, 
and to filter the false positive amplicons.


}
\details{
\tabular{ll}{
Package: \tab AmpliconDuo\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2016-01-14\cr
License: \tab GPL-2 \cr
}

The core of this package is the \code{\link{ampliconduo}} function, that generates for each pair of a split samples an ampliconduo data frame, while statistically analysing the data by Fisher's exact test.
Ampliconduo data frames, or lists of these, are the input required for all other functions of this package.
\itemize{
\item{\code{\link{plotAmpliconduo}}}{ 
plots for an ampliconduo the amplicon frequencies (number of reads per amplicon) of sample A vs. amplicon frequencies of sample B, highlighting amplicons displaying a significant deviation between both samples.
}
\item{\code{\link{plotAmpliconduo.set}}}{
does the same as \code{plotAmpliconduo} but accepts a list of ampliconduo data frames and arranges the plots in a 2-dimensional array.
}
\item{\code{\link{plotORdensity}}}{ 
generates a histogram plot of the amplicon frequency odds ratio density for an ampliconduo data frame. For multiple data frames
organizes the plots in a 2-dimentional array.
}
\item{\code{\link{discordance.delta}}}{
calculates delta (\eqn{\Delta}) and delta prime (\eqn{\Delta'}), the fraction of amplicon frequencies and amplicons, respectively, with a false discovery rate below a certain threshold \eqn{\theta} as a measure of discordance between two amplicon data sets A and B.
}
\item{\code{\link{filter.ampliconduo}}}{
applies filter criteria to an ampliconduo data frame deciding which amplicons are going to be rejected.
}
\item{\code{\link{filter.ampliconduo.set}}}{
same as \code{filter.ampliconduo} for a list af ampliconduo data frames.
}
\item{\code{\link{accepted.amplicons}}}{
returns the indices of those amplicons that have passed the filter criteria.
}
}
}

\author{
Anja Lange (\email{anja.lange@uni-due.de})
and
Daniel Hoffmann (\email{daniel.hoffmann@uni-due.de})

Maintainer: Anja Lange (\email{anja.lange@uni-due.de})

}


\keyword{ package }

\references{
Lange A, Jost S, Heider D, Bock C, Budeus B, Schilling E, Strittmatter A, Boenigk J, Hoffmann D: AmpliconDuo: A Split-Sample Filtering Protocol for High-Throughput Amplicon Sequencing of Microbial Communities. 
\emph{(PLoS One. 2015 Nov 2;10(11))}
}

\examples{

## load test amplicon frequency data ampliconfreqs and vector with sample names site.f
data(ampliconfreqs)
data(site.f)

## generating ampliconduo data frames 
## depending on the size if the data sets, may take some time
ampliconduoset <- ampliconduo(ampliconfreqs[,1:4], sample.names = site.f[1:2])

## plot amplicon read numbers of sample A  vs. amplicon read numbers of sample B,
## indicating amplicons with significant deviations in their occurence across samples
plotAmpliconduo.set(ampliconduoset, nrow = 3)

## calculate discordance between the two data sets of an ampliconduo
discordance <- discordance.delta(ampliconduoset)

## plot the odds ratio density of ampliconduo data
plotORdensity(ampliconduoset)

## apply filter criteria to remove/mark spurious amplicons
ampliconduoset.f <- filter.ampliconduo.set(ampliconduoset, min.freq = 1, q = 0.05)

## return indices of accepted amplicons, indices correspond to indices of the ampliconfreqs data, 
## that were used as input for the ampliconduo function
accep.reads <- accepted.amplicons(ampliconduoset.f)
}
