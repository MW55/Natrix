\name{filter.ampliconduo}
\alias{filter.ampliconduo}

\title{
AmpliconDuo Amplicon Filter
}
\description{

Marks or removes amplicons from an ampliconduo data frame according to the specified filter criteria.

}
\usage{
filter.ampliconduo(x, min.freq = 1, OR = NULL, q = NULL, p = NULL, remove = FALSE)
}
\arguments{
  \item{x}{
 Data frame, an ampliconduo data frame returned by the \code{\link{ampliconduo}} function.
}
  \item{min.freq}{
Optional. Integer, minimium frequency/read count for an amplicon in each of the two amplicon sets to be retained.
Default value is 1.
}
  \item{OR}{
Optional. Numeric, minimum odds ratio for an amplicon to be retained. If no value is specified the odds ratio is 
excluded from the filter criteria.
}
  \item{q}{
Optional. Numeric, minimum value for \code{q}, the adjusted p-value for
an amplicon to pass the filter. If no value 
is specified, \code{q} is excluded from the filter criteria.
}
  \item{p}{
Optional. Numeric, minimum p-value for an amplicon to pass the filter. If no value for \code{p}
is specified \code{p} is excluded from the filter criteria.
}
  \item{remove}{
Optional. Logical, decides whether amplicons that fail the filter criteria should be removed (\code{TRUE}),
or retained (\code{FALSE}).
}
}
\details{

Takes the \code{ampliconduo} data frame \code{x} and applies each filter criterion that is selected to each amplicon.
If an amplicon \emph{i} does not pass each of the applied criteria, the logical value in column \code{rejected} in row \emph{i} is set to \code{TRUE}.
In case the parameter \code{remove} was set to \code{TRUE}, all amplicons with \code{rejected = TRUE} are removed. The position information with respect to the data used as input for the \code{\link{ampliconduo}} call is retained.
}
\value{

Data frame corresponding to the input \code{x}, but with the adjustments in the \code{rejected} column according to the specified filter criteria, or removed rows (\code{removed = TRUE}).
}

\author{
Anja Lange
&
Daniel Hoffmann
}

\seealso{
\code{\link{ampliconduo}}, generates the input data \code{x} for this method.

\code{\link{accepted.amplicons}}, returns the indices of amplicons that
have passed the filter criteria.
}
\examples{
##load example data
data(amplicons)

## extract the first ampliconduo data frame
ampliconduo1 <- amplicons[[1]]

## apply filter criteria 
ampliconduo1.f <- filter.ampliconduo(ampliconduo1)
ampliconduo1.f <- filter.ampliconduo(ampliconduo1, min.freq = 2, remove = TRUE)

## to return a list with the indices (corresponding to the indices of the data 
## the ampliconduo function was called on) of all amplicons that passed the filter criteria 
good.reads <- accepted.amplicons(ampliconduo1.f)

}

\keyword{manip}
