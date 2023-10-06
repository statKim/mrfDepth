\name{sprojmedian}
\alias{sprojmedian}
\title{Location estimates based on skewness-adjusted projection depth}
\description{
Computes a skewness-adjusted projection depth based location estimate of a 
\eqn{p}-dimensional dataset \code{x}.
}
\usage{
sprojmedian(x, sprojection.depths = NULL, options = NULL)
}
\arguments{
    \item{x}{An \eqn{n} by \eqn{p} data matrix with observations
             in the rows and variables in the columns.}
  \item{sprojection.depths}{Vector containing the skewness-adjusted projection
                  depth values of the observations in \code{x}.}
  \item{options}{ A list of options to pass to the \code{sprojdepth} routine.
                  See \code{sprojdepth} for more details. }

}
\details{
 The algorithm depends on the function \code{sprojdepth} to compute the
 skewness-adjusted projection depth values of the observations in \code{x}. 
 If these depth have already been
 computed they can be passed as an optional argument to save computing time.
 If not, the skewness-adjusted projection depth values will be computed and the user 
 can pass a list with options to the \code{sprojdepth} function.

 It is first checked whether the data lie in a subspace of dimension smaller 
 than \eqn{p}. If so, a warning is given, as well as the dimension of the 
 subspace and a direction which is orthogonal to it.
}
\value{
  A list with component: \cr
    \item{max}{The point of \code{x} with maximal skewness-adjusted projection 
               depth. If multiple points have maximum depth, their center
               of gravity is returned.}
}
\seealso{
  \code{\link{adjOutl}}, \code{\link{sprojdepth}}, \code{\link{dirOutl}}, \code{\link{outlyingness}}
}
\references{
Hubert M., Van der Veeken S. (2008). Outlier detection for skewed data. \emph{Journal of Chemometrics}, \bold{22}, 235--246.
 
Hubert M., Rousseeuw P.J., Segaert P. (2015). Multivariate functional outlier detection. \emph{Statistical Methods & Applications}, \bold{24}, 177--202.
}
\author{P. Segaert}
\examples{
# Compute a location estimate of a two-dimensional dataset.
data(bloodfat)

result <- sprojmedian(x = bloodfat)
plot(bloodfat, pch = 16)
points(result$max, col = "red", pch = 18, cex = 1.5)

# Options for the underlying sprojdepth routine may be passed 
# using the options argument. 
result <- sprojmedian(x = bloodfat,
                      options = list(type = "Rotation",
                                     ndir = 1000
                                     )
                     )
plot(bloodfat, pch = 16)
points(result$max, col = "red", pch = 18, cex = 1.5)

# One may also compute the depth values of the observations in the data
# separately. This avoids having to recompute them when computing the median. 
depth.result <- sprojdepth(x = bloodfat)
result <- sprojmedian(x = bloodfat, 
                      sprojection.depths = depth.result$depthX)
}
\keyword{multivariate}
