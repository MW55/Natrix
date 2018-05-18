\name{plotAmpliconduo}
\alias{plotAmpliconduo}
\title{
Frequency-Frequency Plot Focusing On Significant Deviations In Amplicon Occurences
}
\description{
Applied to an ampliconduo data frame, one element of the return value of
the \code{\link{ampliconduo}} function. Generates a plot of \code{freqB}
over \code{freqA} (the read numbers of the same amplicon in both halves
A and B of a split sample). For amplicons that have significantly
deviating read numbers, i.e. with a p-value or adjusted p-value below a
certain treshold, points are colored differently (default: red).
}
\usage{
plotAmpliconduo(x, color.treshold = 0.05, xlab = "Abundance (PCR A)", 
ylab = "Abundance (PCR B)",main = NULL, log = "xy", corrected = TRUE, asp = 1, 
legend.position = NULL, save = FALSE, path = NULL, file.name = NULL,
format = "jpeg", h.start = 0, ...)
}

\arguments{
  \item{x}{
Ampliconduo data frame, an element of the returned list of the \code{\link{ampliconduo}} function.
}
  \item{color.treshold}{
Optional. Numeric value specifying at which p-value or adjusted p-value points in the plot are drawn in complementary color.
Default value is 0.05.
}
  \item{xlab}{
Optional. Character indicating the x-axis label. Default is \dQuote{Abundance (PCR A)}.
}
  \item{ylab}{
Optional. Character indicating the y-axis label. Default is \dQuote{Abundance (PCR B)}.
}
  \item{main}{
Optional. Character specifying the overall title of the plot. If no value is passed, takes the 
sample name from the \code{x$sample} column.
}
  \item{log}{
Optional. Character specifying the variables to transform to log (\dQuote{},\dQuote{x}, \dQuote{y}, or \dQuote{xy}).
Default is \dQuote{xy}.
}
  \item{corrected}{
Optional. Logical to decide whether the p-value (\code{FALSE}), or corrected p-value (\code{TRUE}) is used for coloring.
Default value is \code{TRUE}.
}
  \item{asp}{
Optional. Numeric value, the y/x aspect ratio. Default is 1.
}
  \item{legend.position}{
Optional. Numeric vector of length two. Defines the position of the legend. By default tries to find a position that fits
best the arrangement of the plots.
}
  \item{save}{
Optional. Logical value indicationg if the plot should be saved to file. Default value is \code{FALSE}.
}
  \item{path}{
Optional. Character, in case the argument \code{save} was set to \code{TRUE}, specifies the directory (no backslash or slash at the end) for saving.
By default (parameter \code{save} is set to TRUE) the plot is saved in the working directory.
}
  \item{file.name}{
Optional. If  argument \code{save} was set to \code{TRUE} specifies the file name for the plot.
The default name is \dQuote{ampliconduo_<date>.<file extension>}
}
  \item{format}{
Optional. Character specifying the format of the saved file. One of \dQuote{eps}, \dQuote{ps}, \dQuote{tex}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{tiff}, \dQuote{png}, \dQuote{bmp}, \dQuote{svg} and \dQuote{wmf} (windows only). Default format is \dQuote{jpeg}.
}
  \item{h.start}{
Optional. Numeric value between 0 and 360 defines the color of the plotted points. 
Default value is 0 (blue-green, red). 
}
  \item{\dots}{
Optional. Allows to pass other aesthetics.
}
}

\references{
\url{http://had.co.nz/ggplot2/book/qplot.pdf}
}

\author{
Anja Lange 
&
Daniel Hoffmann
}

\seealso{

\pkg{ggplot2} package

\code{\link[ggplot2]{qplot}} internally used to create the plot.

\code{\link{plotAmpliconduo.set}}, generates a very similar plot for a list of ampliconduo data frames.

\code{\link{ampliconduo}}, generates the input data.
}
\examples{
## load example data
data(amplicons)

## extract the second ampliconduo data frame
ampliconduo2 <- amplicons[[2]]

## plot the amplicon frequencies of the ampliconduo data frame 
plotAmpliconduo(ampliconduo2, main = "ampliconduo_2")
plotAmpliconduo(ampliconduo2, main = "ampliconduo_2", h.start = 50, log = "")
plotAmpliconduo(ampliconduo2, h.start = 50, log = "", asp = 2, corrected = FALSE)
}

\keyword{hplot}
