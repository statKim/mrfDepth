\name{fHeatmap}
\alias{fHeatmap}
\title{Draws a heatmap of functional depth values or distances}
\description{
Draws a heatmap of depth values or distances of functional data.
}
\usage{
fHeatmap(rowValues, cellValues, type, legend.title = "")
}
\arguments{
  \item{rowValues}{Vector of functional depth or distance values. Each value should correspond with the functional depth or distance of an observation from a (multivariate) functional data set.}
  \item{cellValues}{Matrix of multivariate depth or distance values. The value in row \eqn{i} and column \eqn{j} should correspond with the multivariate depth or distance of observation \eqn{i} at time point \eqn{j}.}
  \item{type}{One of \code{"depth"} or \code{"distance"}. Determines whether a depth or a distance heatmap is drawn.}
  \item{legend.title}{Title of the legend.}
}
\details{On the vertical axis the functional data are sorted from top to bottom according to their functional depth or distance value as provided in \code{rowValues}. When \code{type = "depth"}, the deepest observation is put at the bottom. When \code{type = "distance"}, the rows are sorted in decreasing order of their functional distance.

Each cell of the map is colored according to the multivariate depth provided in \code{cellValues}. For a depth-based heatmap, the smallest
depth value is  white and the overall highest depth value is colored dark green. A distance-based heatmap colors the cells from white to dark red.  
}
\references{
Hubert M., Rousseeuw P.J., Segaert P. (2015). Multivariate functional outlier detection (with rejoinder). \emph{Statistical Methods & Applications}, \bold{24}, 177--202.
}

\author{P. Segaert}

\examples{
data(octane)
Result <- mfd(octane, diagnostic = TRUE, type = "sprojdepth")
Plot <- fHeatmap(rowValues = Result$MFDdepthZ,
                 cellValues = Result$crossdepthZ,
                 type = "depth",
                 legend.title = "SPD")
Plot

Result <- fOutl(octane, diagnostic = TRUE, type = "fAO")
Plot <- fHeatmap(rowValues = Result$fOutlyingnessZ,
                 cellValues = Result$crossDistsZ,
                 type = "distance",
                 legend.title = "AO")                 
Plot
}
\keyword{Graphical}
