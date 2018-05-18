\name{discordance.delta}
\alias{discordance.delta}

\title{
Measure Of Discordance Between Two Amplicon Data Sets 
}
\description{
Calculates delta (\eqn{\Delta}) and delta prime (\eqn{\Delta'}), the fraction of amplicon frequencies and amplicons, respectively, with a false discovery rate below a certain threshold \eqn{\theta} as a measure of discordance between two amplicon data sets A and B.

}
\usage{
discordance.delta(x, names = NULL, theta = 0.05, corrected = TRUE,
printToTex = FALSE, directory = NULL, file.name = NULL)
}
\arguments{
  \item{x}{
   A list of amplicon duo data frames as returned by the \code{\link{ampliconduo}} function.
}
  \item{names}{
Optional. Vector or list of characters specifying the sample/amplicon pair names. By default 
names are taken from the element names of \code{x}.
}
  \item{theta}{
Optional. Numeric, threshold for the false discovery rate. Default value is 0.05.
}
  \item{corrected}{
Optional. Logical, indicates whether the p-value from Fisher's exact test (FALSE) or the 
adjusted p-value (TRUE), here called q, is used for calculation of \eqn{\Delta} and \eqn{\Delta'}.
}
  \item{printToTex}{
Optional. Logical, if \code{TRUE} writes a file with the resulting table in Tex format. Default is \code{FALSE}.
}
  \item{directory}{
Optional. If \code{printToTex} is \code{TRUE}, the output file will be saved in the specified directory.
}
  \item{file.name}{
Optional. If \code{printToTex} is \code{TRUE}, the output file will have the name \code{file.name}.
}
}
\details{
Calculates \eqn{\Delta} and \eqn{\Delta'}, the fraction of frequencies of amplicons and amplicons, respectively, with false discovery rate below a certain threshold \eqn{\theta} as a measure of discordance between two amplicon data sets A and B with occurence \eqn{r_{(iAS)}} of amplicon \emph{i} in amplicon set A of sample S (the ampliconduo data frame).
\eqn{\Delta} and \eqn{\Delta'} are defined as follows:

\deqn{ 
  \Delta_{S\theta}  = \frac{\sum_{i=1}^{n_S}(r_{iAS}+r_{iBS})\delta(q_{iS}<\theta)}{\sum_{i=1}^{n_S}(r_{iAS}+r_{iBS})}\\
}{}

\deqn{
  \Delta'_{S\theta} = \frac{\sum_{i=1}^{n_S}\delta(q_{iS}<\theta)}{n_S},\\
}{}

\deqn{
  \mathrm{with\;} \delta(q_{iS}<\theta) =
  \left\{
  \begin{array}{cc}
  1 & \mathrm{\;for\;} q_{iS} < \theta\\
  0 & \mathrm{\;for\;} q_{iS} \ge \theta
  \end{array}
  \right. 
}{}

\if{html}{\figure{discordanceDelta.png}{options: width="400"}}


for number \eqn{n_S} of amplicons detected in sample/ampliconduo S. \eqn{\Delta_{S\theta}} and \eqn{\Delta'_{S\theta}} are located between 0 (no discordance, i.e. no statistically significant deviations between experimental branches) and 1 (complete discordance).

\code{x} is the return value of a \code{ampliconduo} call.

}
\value{
Data frame with three columns. The first column contains the sample/amplicon pair names. Second and third column 
harbor the corresponding \eqn{\Delta} and \eqn{\Delta'} values, respectively.
}

\author{
Anja Lange
&
Daniel Hoffmann
}


\seealso{
 \code{\link{ampliconduo}}, generates the expected data format for \code{x}
 
 \code{\link[xtable]{xtable}}, used to convert the returned data frame into a \code{Latex} table.
}

\examples{
## load example ampliconduo data frame
data(amplicons)

## calculate the discordance between amplicon data sets of an ampliconduo data frame
dd.a <- discordance.delta(amplicons)
dd.b <- discordance.delta(amplicons, theta = 0.1)
dd.c <- discordance.delta(amplicons, printToTex = TRUE)
}


\keyword{htest}
