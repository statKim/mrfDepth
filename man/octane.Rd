\name{octane}
\alias{octane}
\title{Near infrared spectra of gasoline samples}
\description{
Near infrared absorbance spectra (NIR) of 39 gasoline samples over 226 wavelengths ranging from 1102nm to 1552nm with measurements every two nanometers. Six of the samples (25, 26, 36-39) contain added alcohol.
}
\usage{
data(octane)
}
\format{
A three-dimensional \eqn{t = 226} by \eqn{n = 39} by \eqn{p = 1} array.
}
\source{
Esbensen K., Schonkopf S., Midtgaard T. (1994). Multivariate Analysis
in Practice. Trondheim: Camo.
}
\references{
Hubert M., Rousseeuw P.J., Vanden Branden K. (2005).
ROBPCA: a new approach to robust principal component analysis.
\emph{Technometrics}, \bold{47}, 64--79.
}
\examples{
data(octane)
matplot(octane[,1:30,1], type = "l", lty = 1, col = "black")
}
\keyword{datasets}
