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
res = susy(data[, c("var1","var2")], segment=30, fps=15)
names(res)

## compute SUSY for var2-var3
res = susy(data[, c("var2","var3")], segment=30, fps=15)
names(res)

## compute SUSY for var1-var2 and var3-var4
res = susy(data, segment=30, fps=15)
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

## export to flat file
df = as.data.frame(res)
write.table(df, file="correlation.txt", sep=",", row.names=FALSE)
```

[`susy` function manual](https://WTschacher.github.io/SUSY/reference/susy.html)
