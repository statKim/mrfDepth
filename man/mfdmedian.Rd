\name{mfdmedian}
\alias{mfdmedian}
\title{Multivariate functional median for functional data}
\description{
Computes the multivariate functional median, an estimate for the central tendency of multivariate functional data.
}
\usage{
  mfdmedian(x, type = "hdepth", crossDepthsX = NULL,
            depthOptions = NULL, centerOption = "maxdepth")
}
\arguments{
       \item{x}{A three-dimensional \eqn{t} by \eqn{n} by \eqn{p} array,
                with \eqn{t} the number of observed time points,
                \eqn{n} the number of functional observations
                and \eqn{p} the number of measurements
                for every functional observation at every time point.}
    \item{type}{The depth used in the computations.
                One of the following options: \code{"hdepth"},
                \code{"projdepth"}, \code{"sprojdepth"}, \code{"dprojdepth"},
                \code{"sdepth"}. \cr  Defaults to \code{"hdepth"}.}
\item{crossDepthsX}{Depth values at each time point. Can be used to save computing time.}
\item{depthOptions}{A list of options to pass to the function
                    that computes the cross-sectional depths.} 
\item{centerOption}{When equal to \code{"maxdepth"} the functional median equals 
                    at each time point the point with cross-sectional maximal depth. 
                    When \code{type} is equal to \code{"projdepth"}, also a weighted center 
                    of gravity can be computed based on Huber weights (see                                               \code{\link{projmedian}}). Then \code{centerOption} should be set to                                 \code{"Huber"}. Defaults to \code{"maxdepth"}.}
}
\details{
The multivariate functional median of a multivariate functional data set is defined
as the multivariate curve connecting the cross-sectional multivariate depth 
medians (Claeskens et al., 2014). 
The MFD median can be computed in all dimensions \eqn{p} using halfspace depth, 
projection depth, skewness-adjusted projection depth or directional projection depth. 
The simplicial depth can only be used for \eqn{p \le 2}.

It is possible that at certain time points a part of the algorithm can not be 
executed due to e.g. exact fits. In that case the estimate for the center will 
be set to NaN. A warning is issued at the end of the algorithm to signal 
these time points. Furthermore the output contains an extra argument giving the 
indices of the time points where problems occured.
}

\value{
 A  list with the following
  component:
       \item{MFDmedian}{An \eqn{t} by \eqn{p} matrix containing the estimated central curve.}
       \item{IndFlagExactFit}{Vector containing the indices of the
                              time points for which an exact fit is
                              detected.}
}
\references{
Claeskens G., Hubert M., Slaets L., Vakili K. (2014). Multivariate functional halfspace depth.
\emph{Journal of the American Statistical Association,} \bold{109}, 411--423.
}
\seealso{
\code{\link{hdepth}}, \code{\link{projdepth}}, 
\code{\link{sprojdepth}}, \code{\link{dprojdepth}}, \code{\link{sdepth}},
\code{\link{mfd}}
}
\examples{
# Consider a bivariate functional sample.
data(characterA)
Data <- characterA[, 1:50, ]
Result <- mfdmedian(Data)

# Plot the data and the functional median
par(mfrow = c(1,2))
matplot(Data[, , 1], type = "l", col = "black", lty = 1, ylab = "x-coordinate")
matlines(Result$MFDmedian[, 1], type = "l", col = "red", lwd = 2)
matplot(Data[, , 2], type = "l", col = "black", lty = 1, ylab = "y-coordinate")
matlines(Result$MFDmedian[, 2], type = "l", col = "red", lwd = 2)
par(mfrow = c(1,1))

# Other depth functions such as the adjusted outlyingness may also 
# be used to determine the cross-sectional depth median. 
# In this case the depth median is computed by the 
# sprojmedian routine. 
# Optional arguments used by the sprojmedian routine may be specified
# using the depthOptions. For example one might choose the
# "Rotation" option with 300 directions. 
Result <- mfdmedian(Data, type = "sprojdepth",
                    depthOptions = list(type = "Rotation",
                                        ndir = 300))
par(mfrow = c(1,2))
matplot(Data[, , 1], type = "l", col = "black", lty = 1, ylab = "x-coordinate")
matlines(Result$MFDmedian[, 1], type = "l", col = "red", lwd = 2)
matplot(Data[, , 2], type = "l", col = "black", lty = 1, ylab = "y-coordinate")
matlines(Result$MFDmedian[, 2], type = "l", col = "red", lwd = 2)
par(mfrow = c(1,1))

# If the user already placed a call to the mfd routine
# with the diagnostic options set to TRUE, the 
# mfdmedian can easily be obtained by passing the cross-sectional 
# depths. This considerably saves computing time.  
tResult <- mfd(x = Data, type = "sprojdepth", diagnostic = TRUE)
Result <- mfdmedian(Data, type = "sprojdepth",
                    crossDepthsX = tResult$crossdepthX,
                    )
  
# Univariate curves should also be represented as arrays
Data.x <- Data[, , 1, drop = FALSE]
Result <- mfdmedian(Data.x)
matplot(Data.x[ , , 1], type = "l", col = "black", lty = 1, ylab = "x-coordinate")
matlines(Result$MFDmedian[, 1], type = "l", col = "red", lwd = 2)
}
\author{P. Segaert, M. Hubert}
\keyword{functional}
