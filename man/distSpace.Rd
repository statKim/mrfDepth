\name{distSpace}
\alias{distSpace}
\title{distSpace}
\description{
Computation of distance space representation.
}
\usage{
distSpace(trainingData, testData = NULL, type = "bagdistance", options = NULL)
}
\arguments{
  \item{trainingData}{A list of \eqn{n} by \eqn{p} matrices containing  multivariate data or a list of \eqn{t} by \eqn{n} by \eqn{p} arrays containing functional data. Each element of the list contains the observations from one group. The training data should contain at least two groups.}
  \item{testData}{An \eqn{m} by \eqn{p} matrix containing multivariate test data or a \eqn{t} by \eqn{m} by \eqn{p} array for functional data.}
  \item{type}{The distance used in the computations.
                For multivariate data one of the following options: \code{bagdistance},
                \code{outlyingness}, \code{adjOutl} or \code{dirOutl}. 
                Defaults to \code{bagdistance}.\cr
                For functional data one of the following options: \code{fbd}, \code{fSDO}, 
                \code{fAO} or \code{fDO}. Defaults to \code{fbd}.
                }
  \item{options}{A list of options to pass to the function
                  computing the underlying distance. \cr
                  See \code{bagdistance}, \code{outlyingness}, \code{adjOutl},
                  \code{dirOutl} or \code{fOutl} for more information.}
}
\details{
The distance space representation is a tool in supervised classification and was introduced in Hubert et al. (2016) as a generalisation of the depth-depth representation of a multivariate sample. Based on a distance transform, an observation (be it multivariate or functional) is mapped to its representation in distance space. The distance transformation consists of mapping the observation to a vector containing at coordinate \eqn{i} the distance to the training group \eqn{i}. After transformation, any multivariate classifier may be used to classify new observations in distance space. Typically the \eqn{k}-nearest neighbour algorithm is used. 

Different options are available to compute the distance to each of the training groups. For multivariate data, the user may choose between the bagdistance or any of the projection type distances including the Stahel-Donoho outlyingness, the adjusted outlyingness or the directional outlyingness. For functional data, the user may opt to employ the functional bagdistance (fbd), the functional Stahel-Donoho outlyingness (fSDO), the functional skweness-adjusted outlyingness (fAO) or the functional directional outlyingness (fDO). Options available in each of the underlying distance routines may be passed down using the \code{options} argument.
}
\value{
A \eqn{q} by \eqn{(p+1)} matrix composed of two blocks. The first block contains the observations in the training set (rows) with in each column the distance to each of the groups. The last column contains a label indicating the original group membership of the observation. The second block contains the observations in the test set, if any, with in each column the distance to the different training groups. The last column contains an indicator signaling the observation was part of the test set. 
}
\references{
Hubert M., Rousseeuw P.J., Segaert P. (2017). Multivariate and functional classification using depth and distance. \emph{Advances in Data Analysis and Classification}, \bold{11}, 445--466.
}
\author{P. Segaert}

\examples{
data(plane)

# Build the training data
Mirage <- plane$plane1[, 1:25, 1, drop = FALSE]
Eurofighter <- plane$plane3[, 1:25, 1, drop = FALSE]
trainingData <- list(group1 = Mirage,
                     group2 = Eurofighter)

# Build the test data 
Mirage.t <- plane$plane1[, 26:30, 1, drop = FALSE]
Eurofighter.t <- plane$plane3[, 26:30, 1, drop = FALSE]
testData <- abind::abind(Mirage.t, Eurofighter.t, along = 2) 

# Transform the data into distSpace
Result <- distSpace(trainingData = trainingData, testData = testData, type="fbd")

# Plot the results
plotColors <- c(rep("orange", dim(Mirage)[2]),
                rep("blue", dim(Eurofighter)[2]),
                rep("green3", dim(testData)[2]))
plot(Result[, 1:2, ],
     col = plotColors, pch=16,
     xlab = "distance to Mirage", ylab = "distance to Eurofighter",
     main = "distSpace representation of Mirage and Eurofighter")
     legend("bottomleft", legend = c("Mirage","Eurofighter", "test data"), pch = 16, 
       col = c("orange","blue", "green3"))
       
}

\keyword{functional}
\keyword{multivariate}
