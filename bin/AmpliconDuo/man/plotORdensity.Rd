\name{plotORdensity}
\alias{plotORdensity}

\title{
Probability Density Plot Of Odds Ratios
}
\description{
Plots for an ampliconduo data frame probability densities of the odds ratios of amplicon occurences in the two amplicon data sets.
The function allows to shift the two extrema (odds ratios OR = 0 and OR
= infinity) to the edges of the plot.
Plots of multipe ampliconduo data frames are arranged in a 2-dimensional array with shared scales.

}
\usage{
plotORdensity(x, log = "x", ncol = 2, adjust.zeroinf = TRUE, zero.pos = 0.005, 
inf.pos = 200, binwidth = 0.15, color = "black", xlab = "odds ratio", 
save = FALSE, path = NULL, file.name = NULL, format = "jpeg", ...)
}

\arguments{
  \item{x}{
List or a single ampliconduo data frame, return value of the \code{\link{ampliconduo}} function.
}
  \item{log}{
Optional. Character specifying the variables to transform to log (\dQuote{},\dQuote{x}, \dQuote{y}, or \dQuote{xy}).
Default is \dQuote{x}.
}
  \item{ncol}{
Optional. Integer value specifying the numer of columns used to arrange the plots. Default is 2.
}
  \item{adjust.zeroinf}{
Optional. Logical, specifies whether the density bar for 0 and inf should be shifted. Default value is \code{TRUE}.
}
  \item{zero.pos}{
Optional. Numeric, in case \code{adjust.zeroinf} is set to \code{TRUE}, determines the position for the odds ratio = 0 bar.
Default value is 0.05.
}
  \item{inf.pos}{
Optional. Numeric, in case \code{adjust.zeroinf} is set to \code{TRUE}, gives the position for the odds ratio = inf bar.
Default value is 200.
}
  \item{binwidth}{
Optional. Numeric, bin width to use, default is 0.15.
}
  \item{color}{
Optional. Character, name of the color used to draw the density bars. Default is \dQuote{black}.
}
  \item{xlab}{
Optional. Character, label for the x-axis. Default is \dQuote{odds ratio}.
}
  \item{save}{
Optional. Logical, \code{TRUE} to save the plot, else \code{FALSE}. Default is \code{FALSE}.
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
  \item{\dots}{
Optional. Allows to pass other aesthetics.
}
}


\author{
Anja Lange 
&
Daniel Hoffmann
}

\seealso{
\code{\link[ggplot2]{qplot}}, used by \code{\link{plotAmpliconduo.set}} to create the plots.

\code{\link[ggplot2]{facet_wrap}}, called for 2-dimensional arrangement of the plots.

\code{\link{ampliconduo}}, generates the input data.
}

\examples{
## loads example data of ampliconduo data frames
data(amplicons)

## plot odds ratio density for amplicon frequencies in ampliconduo data frames
plotORdensity(amplicons)
plotORdensity(amplicons[1:4], binwidth = 0.1, color = "magenta")
plotORdensity(amplicons[[1]], binwidth = 0.1, color = "orange", main = "Sample FU25")
plotORdensity(amplicons[1:4],  color = "darkblue", ncol = 2, save = TRUE)

}

\keyword{hplot}

