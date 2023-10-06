\name{projmedian}
\alias{projmedian}
\title{Location estimates based on projection depth}
\description{
Computes a projection depth based location estimate of a \eqn{p}-dimensional 
dataset \code{x}.
}
\usage{
projmedian(x, projection.depths = NULL, options = NULL)
}
\arguments{
    \item{x}{An \eqn{n} by \eqn{p} data matrix with observations
             in the rows and variables in the columns.}
  \item{projection.depths}{Vector containing the projection depth of the observations in \code{x}.}
  \item{options}{A list of options to pass to the \code{projdepth} routine.
                 See \code{projdepth} for more details. }

}
\details{
 The algorithm depends on the function \code{projdepth} to compute the
 projection depth of the observations in \code{x}. If these depth values have already been computed, they can be passed as an optional
 argument to save computing time. If not, the projection depth values will be
 computed and the user can pass a list with options to the
 \code{projdepth} function.

 It is first checked whether the data lie in a subspace of dimension smaller 
 than \eqn{p}. If so, a warning is given, as well as the dimension of the 
 subspace and a direction which is orthogonal to it.

}
\value{
  A list with components: \cr
     \item{max}{The point of \code{x} with maximal projection depth.
                If multiple points have maximum depth, their center
                of gravity is returned.}
     \item{Huber}{A weighted center of gravity of all observations.
                  The weights are defined by the Huber
                  function with parameter
                        \code{alpha = 1 / (1+sqrt(qchisq(0.95, p)))}
                  Observations for which the projection depth
                  is more than alpha receive weight 1,
                  other points receive weight
                  \code{( sqrt(qchisq(0.95, p)) / outlyingness )^2.}}
}
\seealso{
  \code{\link{outlyingness}}, \code{\link{projdepth}}, \code{\link{adjOutl}}, \code{\link{dirOutl}}
}
\references{
Maronna, R.A., Yohai, V.J. (1995). The behavior of the Stahel-Donoho robust multivariate estimator. \emph{Journal of the American Statistical Association}, \bold{90}, 330--341.
}
\author{P. Segaert}
\examples{
# Compute a location estimate of a two-dimensional dataset.
if (requireNamespace("robustbase", quietly = TRUE)) {
  BivData <- log(robustbase::Animals2)
} else {
  BivData <- matrix(rnorm(120), ncol = 2)
  BivData <- rbind(BivData, matrix(c(6, 6, 6, -2), ncol = 2))
}

result <- projmedian(x = BivData)
plot(BivData, pch = 16)
points(result$max, col = "red", pch = 18, cex = 1.5)
points(result$Huber, col = "blue", pch = 3, cex = 1.5)

# Options for the underlying projdepth routine may be passed 
# using the options argument. 
result <- projmedian(x = BivData,
                     options = list(type = "Rotation",
                                    ndir = 100,
                                    stand = "unimcd",
                                    h = 0.75*nrow(BivData)))
              
plot(BivData, pch = 16)
points(result$max, col = "red", pch = 18, cex = 1.5)
points(result$Huber, col = "blue", pch = 3, cex = 1.5)

# One may also compute the depth values of the observations in the data
# separately. This avoids having to recompute them when computing the median. 
depth.result <- projdepth(x = BivData)
result <- projmedian(x = BivData, 
                     projection.depths = depth.result$depthX)
}
\keyword{multivariate}
