\name{dprojmedian}
\alias{dprojmedian}
\title{Location estimates based on directional projection depth}
\description{
Computes a directional projection depth based location estimate of a 
\eqn{p}-dimensional dataset \code{x}.
}
\usage{
dprojmedian(x, dprojection.depths = NULL, options = NULL)
}
\arguments{
    \item{x}{An \eqn{n} by \eqn{p} data matrix with observations
             in the rows and variables in the columns.}
  \item{dprojection.depths}{Vector containing the directional projection
                  depth of the observations in \code{x}.}
  \item{options}{A list of options to pass to the \code{dprojdepth} routine.
                  See \code{dprojdepth} for more details. }

}
\details{
 The algorithm depends on the function \code{dprojdepth} to compute the
 directional projection depth of the observations in \code{x}. If these depth values  have already been computed they can be passed as an optional argument to save computing time. If not, directional projection depth values will be computed and the user 
 can pass a list with options to the \code{dprojdepth} function.

 It is first checked whether the data lie in a subspace of dimension smaller 
 than \eqn{p}. If so, a warning is given, as well as the dimension of the 
 subspace and a direction which is orthogonal to it.
}
\value{
  A list with component: \cr
    \item{max}{The observation of \code{x} with maximal directional projection 
               depth. If multiple points have maximum depth, their center
               of gravity is returned.}
   }
\seealso{
  \code{\link{dirOutl}}, \code{\link{dprojdepth}}, \code{\link{adjOutl}}, \code{\link{outlyingness}}
}
\references{
Rousseeuw P.J., Raymaekers J., Hubert M. (2018). A measure of directional outlyingness with applications to image data and video. \emph{Journal of Computational and Graphical Statistics}, \bold{27}, 345--359.
}
\author{J. Raymaekers}
\examples{
# Compute a location estimate of a two-dimensional dataset.
data(bloodfat)
result <- dprojmedian(x = bloodfat)
plot(bloodfat, pch = 16)
points(result$max, col = "red", pch = 18, cex = 1.5)

# Options for the underlying dprojdepth routine may be passed 
# using the options argument. 
result <- dprojmedian(x = bloodfat, options = list(type = "Rotation", ndir = 100))
plot(bloodfat, pch = 16)
points(result$max, col = "red", pch = 18, cex = 1.5)

# One may also compute the depth values of the observations in the data
# separately. This avoids having to recompute them when conmputing the median. 
depth.result <- dprojdepth(x = bloodfat)
result <- dprojmedian(x = bloodfat, dprojection.depths = depth.result$depthX)
}
\keyword{multivariate}
