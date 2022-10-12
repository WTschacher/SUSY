\name{print.susy}
\alias{print.susy}
\title{ susy print method }
\description{
  Prints information about an \code{susy} object.
}
\usage{
  \method{print}{susy}(x, corr.without.amount=TRUE, \dots)
}
\arguments{
  \item{x}{ A \code{susy} object. }
  \item{corr.without.amount}{ Logical, defaults to \code{TRUE} display correlation without the amount. }
  \item{\dots}{ Ignored. }
}
\value{
  Returns \code{x} invisibly.
}
\seealso{\code{\link{susy}}}
\examples{
n = 1000
data = data.frame(
  var1 = runif(n,300,330),
  var2 = runif(n,300,330)
)
res = susy(data)
res
}
\keyword{ data }
