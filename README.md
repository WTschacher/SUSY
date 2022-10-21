SUSY
----

SUSY computes synchrony as windowed cross-correlation based on two-dimensional time series, as described in Tschacher, Rees & Ramseyer (2014).

[R package website](https://WTschacher.github.io/SUSY)

----

Installation
----

```r
install.packages("SUSY", repos="https://WTschacher.github.io/SUSY")
```

Usage
----


```r
library(SUSY)

## random data
n = 10000
data = data.frame(
 var1 = runif(n, 300, 330),
 var2 = runif(n, 300, 330),
 var3 = runif(n, 300, 330),
 var4 = runif(n, 300, 330)
)

## compute SUSY for var1-var2
res = susy(data[, c("var1","var2")], segment=30L, fps=15L)
names(res)

## compute SUSY for var2-var3
res = susy(data[, c("var2","var3")], segment=30L, fps=15L)
names(res)

## compute SUSY for var1-var2 and var3-var4
res = susy(data, segment=30L, fps=15L)
names(res)

## print all SUSY computations
res

## subset (and print) susy object to single results
res["var1-var2"]
res[1]

## plot all SUSY computations, plot type 1
plot(res, type=1)

## plot only first SUSY computations, plot type 1 and 4
plot(res[1], type=c(1,4))

## plot only second SUSY computations, plot type 1, 2, 3, 4, 5
plot(res[2], type=1:5)

## compute SUSY for all permutations of var1, var2, var3 and var4
res = susy(data, permutation=TRUE)
names(res)

## print legacy style
print(res, legacy=TRUE)

## export
df = data.frame(
  xdata = seq(from=-res[[1L]]$params$maxlag , to=res[[1L]]$params$maxlag),
  meanccorrPseudo = res[[1L]]$lagtimes2.data$meanccorrPseudo, meanccorrReal = res[[1L]]$lagtimes2.data$meanccorrReal
)
write.table(df, file="CrossCorrelations.txt", sep=",", row.names=FALSE, col.names=TRUE)
## cleanup
file.remove("CrossCorrelations.txt")
```

[`susy` function manual](https://WTschacher.github.io/SUSY/reference/susy.html)
