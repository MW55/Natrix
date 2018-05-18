\name{accepted.amplicons}
\alias{accepted.amplicons}

\title{ 
List Of Amplicon Indices That Passed Filter Criteria
}
\description{
Returns the indices of those amplicons in an ampliconduo data frame, that passed the applied filter criteria (ampliconduo data frames are filtered using the \code{\link{filter.ampliconduo}} or \code{\link{filter.ampliconduo.set}} function). 
}
\usage{
accepted.amplicons(x)
}

\arguments{
  \item{x}{
An ampliconduo data frame or a list of ampliconduo data frames.
}
}
\details{
Calling this function on an ampliconduo data frame, or a list of the latter, returns the indices of amplicons that passed the applied filter criteria. For each ampliconduo data frame an integer vector is created, and if a list of ampliconduo data frames is supplied with \code{x}, these are pooled in a list.
The returned indices correspond to the data originally used to generate the ampliconduo data frames (parameter \code{A} and \code{B} in  the \code{\link{ampliconduo}} function call).
}
\value{
If \code{x} is an ampliconduo data frame, an integer vector is returned.

In case \code{x} is a list of ampliconduo data frames, a list of integer
vectors is returned, one for each data frame.
}

\author{
Anja Lange
&
Daniel Hoffmann
}

\seealso{
\code{\link{filter.ampliconduo}} and \code{\link{filter.ampliconduo.set}}
}

\examples{
## load example data
data(amplicons)

## apply filter criteria
ampliconduos.f <- filter.ampliconduo.set(amplicons, q = 0.05)

## return a list with accepted amplicons
good.reads <- accepted.amplicons(ampliconduos.f)
}

\keyword{ list }

