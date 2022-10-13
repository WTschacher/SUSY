\name{print.susy}
\alias{print.susy}
\title{ susy print method }
\description{
  Prints information about an \code{susy} object.
}
\usage{
  \method{print}{susy}(x, corr.without.amount=TRUE, legacy=FALSE, \dots)
}
\arguments{
  \item{x}{ A \code{susy} object. }
  \item{corr.without.amount}{ Logical, defaults to \code{TRUE} display correlation without the amount. }
  \item{legacy}{ Logical, defaults to \code{FALSE}, when \code{TRUE} print will produce an output that matches legacy SUSY computation \link{\code{old.susy}}. }
  \item{\dots}{ Extra arguments passed to \code{print.data.frame} method. }
}
\value{
  Returns \code{x} invisibly. Display output to console as a side effect.
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
print(res, corr.without.amount=FALSE)
print(res, digits=4)
print(res, legacy=TRUE)
}
\keyword{ data }
