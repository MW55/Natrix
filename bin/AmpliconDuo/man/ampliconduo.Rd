\name{ampliconduo}
\alias{ampliconduo}
\title{
  Apply Fisher's Exact Tests To Two Amplicon Frequency Sets Of The Same Sample
}
\description{

Implements Fisher's exact test to detect amplicons with significant
deviating read numbers between two amplicon sets of the same sample. The
p-values of the Fisher's exact test are corrected for multiple testing
by computation of the false discovery rates q. This function is intended
to help identifying reads that may be the results of experimental artefacts.
(The calculation can take some time depending on the size of the data sets and 
the computing power.)

}

\usage{
ampliconduo(A, B = NULL, sample.names = NULL, correction = "fdr", ...)
}

\arguments{
  \item{A}{
  A list or a data frame containing amplicon occurences / number of reads per amplicon (integer values). 
}
  \item{B}{
 Optional. A list or a data frame containing amplicon occurences. 
}
  \item{sample.names}{
 Optional. A vector or list of characters with names for the amplicon pairs.
}
  \item{correction}{
 Optional. Specifies the correction method for the p-values from Fisher's exact test.
 Accepts one of the following characters: \code{"holm", "hochberg", "hommel",}
 \code{"bonferroni", "BH", "BY","fdr"} and \code{"none"}. Default is \code{"fdr"}.
 For more details see \code{\link[stats]{p.adjust}}.

}
  \item{\dots}{
Arguments passed to the internally called \code{\link[stats]{fisher.test}} function.
}
}
\details{
If only \code{A} is specified, it is assumed that the list elements 1 &
  2, 3 & 4 etc. of \code{A} are amplicon data of the same sample. In case \code{A} and \code{B}
  are specified, the \emph{i}th frequency set of \code{A} and \code{B}
  are combined. For each amplicon data pair, frequencies at the corresponding
  positions in the lists are assumed to belong to the same amplicon. It is required, that two frequency sets that belong to the same sample, an ampliconduo, have the same length. The \code{ampliconduo} function iterates over all amplicon pairs and performs the following tasks:
\itemize{
\item
  {amplicons with frequency zero in both samples are removed. Position information is retained.}
\item
  {For each amplicon Fisher's exact test using the method \code{\link[stats]{fisher.test}}
  is performed. The p-value, odds ratio and confidence interval are returned. Via the \code{\dots}, arguments
  \code{conf.level}, \code{or} and \code{alternative} can be passed to the \code{\link{fisher.test}} function call. Default values are \code{conf.level} = 0.95, \code{or} = 1 and \code{alternative} = "two.sided".}
\item
  {The p-values are corrected using the \code{\link{p.adjust}} function. By default the method by Benjamini & Hochberg (1995) is used.
  Setting the \code{correction} argument to any of the following characters \code{"holm", "hochberg", "hommel", "bonferroni", "BH", "BY",    "fdr", "none"}, the adjustment method for the p-values can be changed. See function \code{\link{p.adjust}}.
  }
  }

The \pkg{AmpliconDuo} package implements further methods to visualize and filter the returned ampliconduo data frames.
}

\value{
A list of data frames, one for each amplicon pair, that will be called \dfn{ampliconduo data frame} in the following. List entries are named according to the specified \code{sample.names} or numbered.

Each ampliconduo data frame has 9 columns
\itemize{
  \item
  freqA: frequencies of amplicon set A
  \item
  freqB: frequencies of amplicon set B (taken from argument \code{B} if specified)
  \item
  p: p-values calculated with Fisher's exact test
  \item
  OR: odds ratio calculated with Fisher's exact test
  \item
  CI.low: lower confidence limit for \code{OR}
  \item
  CI.up: upper confidence limit for \code{OR}
  \item
  rejected: logical, indicating whether the amplicon was rejected
  \item
  sample: sample name taken from \code{sample.name} if specified, same for all rows in a given data frame 
}

}
\references{
Y Benjamini and Y Hochberg. Controlling the False Discovery Rate: A
Practical and Powerful Approach to Multiple Testing. \emph{Journal of the Royal
Statistical Society}. Series B (Methodological), \strong{57(1)}:289-300, 1995.

}
\author{
Anja Lange
and
Daniel Hoffmann
}

\seealso{
\code{\link[stats]{fisher.test}}, used to calculate the p-value, odds ratio and confidence interval; 

\code{\link[stats]{p.adjust}}, called to correct the p-values;

methods to visualize or further manipulate the ampliconduo data frames:

\code{\link{plotAmpliconduo.set}}, 
\code{\link{plotAmpliconduo}}, 
\code{\link{discordance.delta}}, 

}
\examples{

## loads read numbers from example amplicon data sets
data(ampliconfreqs)
data(site.f)

## generate ampliconduo data frames 
ampliconduos.a <- ampliconduo(A = ampliconfreqs[,1:4], sample.names = site.f[1:2])
ampliconduos.b <- ampliconduo(A = ampliconfreqs[c(1,3)],
B = ampliconfreqs[c(2,4)], sample.names = site.f[1:2],
conf.level = 0.9)

## frequency plot
plotAmpliconduo.set(ampliconduos.a)
}

\keyword{htest}
