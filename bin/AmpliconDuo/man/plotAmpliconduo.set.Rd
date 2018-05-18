\name{plotAmpliconduo.set}
\alias{plotAmpliconduo.set}
\title{
Frequency-Frequency Plot Arrangement Focusing On Significant Deviations In Amplicon Occurences
}
\description{
Called on the return value of the \code{\link{ampliconduo}} function, a list of ampliconduo data frames. Generates for each ampliconduo data frame a plot with \code{freqB} over \code{freqA} and arranges them in a 2-dimensional array, whereas plots in the same row and column share the same scale. Points with a p-value or adjusted p-value below a certain treshold are colored differently (default: red) indicating significant deviations of amplicon occurences between the two samples in an ampliconduo data frame.
}

\usage{
plotAmpliconduo.set(x, color.treshold = 0.05, xlab = "Abundance (PCR A)", 
ylab = "Abundance (PCR B)",log = "xy", corrected = TRUE, asp = 1, nrow = 1,
legend.position = NULL, save = FALSE, path = NULL, file.name = NULL,
format = "jpeg", h.start = 0, ...)
}

\arguments{
  \item{x}{
List of ampliconduo data frames, return value of the \code{\link{ampliconduo}} function.
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
  \item{nrow}{
Optional. Integer value specifying the numer of rows used to arrange the plots. Default is 1.
}
  \item{legend.position}{
Optional. Numeric vector of length two. Defines the position of the legend. By default tries to find a position that fits
the arrangement of the plots best.
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
Optional. Numeric value between 0 and 360, defines the color of the plotted points. Default value is 0 (blue-green, red). 
}
  \item{\dots}{
Optional. Allows to pass other aesthetics.
}

}
\details{
Generates an arrangement of plots from the return value of the \code{\link{ampliconduo}} function, that nicely visualizes those amplicons with a significant deviations in read numbers between the two amplicon data sets. The data in \code{x} are transformed and passed to the \code{\link[ggplot2]{qplot}} function. The 2-dimensional arrangement of the different plots is achieved using \code{\link[ggplot2]{facet_wrap}}. Important aestetic parameters like color, aspect ratio, legend position \dots are easily customized. Optionally, the plot can be saved in a variety of formats.
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

\code{\link[ggplot2]{qplot}}, used by \code{\link{plotAmpliconduo.set}} to create the plots.

\code{\link[ggplot2]{facet_wrap}}, called for 2-dimensional arrangement of the plots.

\code{\link{plotAmpliconduo}}, generates a very similar plot for a single ampliconduo data frame.

\code{\link{ampliconduo}}, generates the input data, an ampliconduo data frame.

}
\examples{
## loads example data of ampliconduo data frames
data(amplicons)

## plot amplicon frequencies of multiple ampliconduo data frames
plotAmpliconduo.set(amplicons[1:4], nrow = 3, h.start = 100)

plotAmpliconduo.set(amplicons[1:4], nrow = 1, corrected = FALSE, color.treshold = 0.1)
}

\keyword{hplot}