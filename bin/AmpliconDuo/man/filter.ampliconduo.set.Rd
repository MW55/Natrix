\name{filter.ampliconduo.set}
\alias{filter.ampliconduo.set}

\title{
AmpliconDuo List Amplicon Filter
}

\description{
Marks or removes amplicons from each ampliconduo data frame in a list according to the specified filter criteria.
}

\usage{
filter.ampliconduo.set(x, min.freq = 1, OR = NULL, q = NULL, p = NULL, remove = FALSE)
}

\arguments{
  \item{x}{
List of ampliconduo data frames, return value of an \code{\link{ampliconduo}} call.
}
  \item{min.freq}{
Optional. Integer, minimium frequency/read count for a given amplicon in each of the two amplicon sets of an ampliconduo to be retained.
Default value is 1.
}
  \item{OR}{
Optional. Numeric, minimum odds ratio for an amplicon to be retained. 
If no value is specified the odds ratio is excluded from the filter criteria.
}
  \item{q}{
Optional. Numeric, minimum value for \code{q}, the adjusted p-value foran amplicon to pass the filter. 
If no value is specified, \code{q} is excluded from the filter criteria.
}
  \item{p}{
Optional. Numeric, minimum p-value for an amplicon to pass the filter.
If no value for \code{p} is specified \code{p} is excluded from the filter criteria.
}
  \item{remove}{
Optional. Logical, decides whether amplicons that fail the filter criteria should be removed (\code{TRUE}),
or retained (\code{FALSE}).
}
}
\details{
For every \code{ampliconduo} data frame in argument \code{x}, applies
  each filter criterion that was specified to each amplicon.
If an amplicon \emph{i} fails any of the applied criteria, the logical value in column \code{rejected} in row \emph{i} is set to \code{TRUE}.
In case the parameter \code{remove} was set to \code{TRUE}, all amplicons with \code{rejected = TRUE} are removed. The position information in respect to the data used as input for the \code{\link{ampliconduo}} call are kept.
This method uses the function \code{\link{filter.ampliconduo}}.
}

\value{
List of ampliconduo data frames. Same as input parameter \code{x} but with the adjustments in the \code{rejected} column according to the specified filter criteria, or removed rows (\code{removed = TRUE})
}

\author{
Anja Lange
&
Daniel Hoffmann
}

\seealso{
\code{\link{filter.ampliconduo}}, performs filtering on single ampliconduo data.frames, is called by this method.

\code{\link{ampliconduo}}, generates the input data \code{x} for this method.

\code{\link{accepted.amplicons}}, returns the indices of amplicons that
have passed the filter criteria.
}
\examples{
## load example data
data(amplicons)

## apply filter criteria
ampliconduos.f <- filter.ampliconduo.set(amplicons)
ampliconduos.f <- filter.ampliconduo.set(amplicons, min.freq = 3, remove = TRUE)

## to return a list with the indices (corresponding to the indices of the data 
## the ampliconduo function was called on) of all amplicons that passed the filter criteria 
good.reads <- accepted.amplicons(ampliconduos.f)
}

\keyword{manip}
