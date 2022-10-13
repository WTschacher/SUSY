\name{susy}
\alias{susy}
\title{ Surrogate Synchrony }
\description{
  Cross-correlations are computed up to a specific lag in seconds \code{Maxlag}, then aggregated within a chosen \code{Segment} of e.g. 30s. Aggregation is performed by transforming correlations to Fisher's Z, then computing mean Z in each segment, then across all segments of the time series. Segment shuffling is used to create surrogate time series, on which the same computations are run. This provides effect sizes \code{ES}. SUSY provides two different synchrony measures of each twin time series: \code{mean Z} and \code{ES of mean Z}; \code{mean absolute_Z} and \code{ES of mean absolute_Z}.\cr
  Legacy implementation is available at \code{\link{old.susy}}.
}
\usage{
susy(x, permutation=FALSE, epoch=30L, fps=15L, maxlag=3L*fps,
  pseudo.simplify=FALSE, pseudo.total=500L)
}
\arguments{
  \item{x}{ A data.frame of numeric columns. }
  \item{permutation}{ Logical, default \code{FALSE} requires \code{x} to have even number of columns which are processed in pairs (1-2, 3-4, etc.). When \code{permutation} is \code{TRUE} then function computes all pairs combinations between columns provided in \code{x} (\code{n*(n-1)/2} pairs). }
  \item{epoch}{ Integer, segment size in seconds, default \code{30L}. Must not be smaller than time window (\code{fps}), also \code{epoch*fps} must not be greater than the timepoints (\code{nrow(x)/2}). Example 1: timepoints: 1000, fps = 10, so the epoch must not be greater than 50s (\code{=1000/(10*2)}); Example 2: maxlag : 90, fps = 10, so epoch must not be less than 9s (\code{=90/10}). }
  \item{fps}{ Integer, frames per second (sampling rate), default \code{15L}. }
  \item{maxlag}{ Integer, window size for \code{\link[stats]{ccf}}, maximum lag in seconds, default \code{3L*fps}. }
  \item{pseudo.simplify}{ Logical, default \code{FALSE}. Simplify calculation of the pseudo data. }
  \item{pseudo.total}{ Integer, default \code{500L}. Ignored when \code{pseudo.simplify} is \code{FALSE} (default). }
}
\details{
  Segments are non-overlapping, and the number of segments that fit into the time series may have a remainder (usually a few seconds at the end of the time series), which is not considered.
}
\value{
  Object of class \code{susy} is returned. Each cross correlation pair is an element in resulting object.
}
\seealso{\code{\link{plot.susy}}, \code{\link{print.susy}}, \code{\link{old.susy}}}
\examples{
n = 1000
data = data.frame(
  var1 = runif(n,300,330),
  var2 = runif(n,300,330),
  var3 = runif(n,300,330)
)

## use only first two columns
res = susy(data[, 1:2])
length(res)
names(res)

## use all columns and permutation
res = susy(data, permutation=TRUE)
length(res)
names(res)

## print susy
res
print(res, legacy=TRUE)

## plot susy
plot(res)
plot(res, type=1:2)
}
\keyword{ data }
