source('fast_sampling.R')

# filter a set of constraints
filterConstraints <- function(constr, sel) {
  list(constr = constr[['constr']][sel, , drop=FALSE],
       rhs = constr[['rhs']][sel],
       dir = constr[['dir']][sel])
}

## Simple 3-dimensional example

n <- 3
constr <- mergeConstraints(
  simplexConstraints(n=n),
  ordinalConstraint(n=n, i=1, j=2),
  ordinalConstraint(n=n, i=1, j=3))

transform <- simplex.createTransform(n=n)
transform.inv <- simplex.createTransform(n=n, inverse=TRUE)

simplex <- cbind(diag(3), rep(1,3)) %*% t(transform.inv)
plot(simplex[,1], simplex[,2])
lines(simplex[c(1,2,3,1),1],simplex[c(1,2,3,1),2], col="black")

constr <- transformConstraints(transform, filterConstraints(constr, constr$dir == "<="))
for (i in 1:nrow(constr$constr)) {
  a <- constr$constr[i,1:2]
  b <- constr$rhs[i] - constr$constr[i,3]
  if (a[1] == 0) {
    abline(h=b/a[2])
  } else {
    abline(a=b/a[2], b=-a[1]/a[2])
  }
}

bb <- createBoundBox(constr, homogeneous=TRUE)
samples.rej <- bbReject(bb$lb, bb$ub, constr, 1000, homogeneous=TRUE)$samples
points(samples.rej[,1],samples.rej[,2],pch="+",col="red")

samples.sds <- sds(constr, 1000, homogeneous=TRUE)$samples
points(samples.sds[,1],samples.sds[,2],pch="+",col="blue")

library(uniformity)
cat(paste0("P-value for uniformity: ", pnorm(testUniformMST(bench=samples.rej, test=samples.sds)), "\n"))

## Six-dimensional example

n <- 6
constr <- mergeConstraints(
  simplexConstraints(n=n),
  ordinalConstraint(n=n, i=2, j=1),
  ordinalConstraint(n=n, i=3, j=1),
  ordinalConstraint(n=n, i=4, j=1),
  ordinalConstraint(n=n, i=5, j=1),
  ordinalConstraint(n=n, i=6, j=1))

transform <- simplex.createTransform(n=n)
constr <- transformConstraints(transform, filterConstraints(constr, constr$dir == "<="))

bb <- createBoundBox(constr, homogeneous=TRUE)
print(system.time(
samples.rej <- bbReject(bb$lb, bb$ub, constr, 1E4, homogeneous=TRUE)$samples
))

print(system.time(
samples.sds <- sds(constr, 1E4, homogeneous=TRUE)$samples
))

library(uniformity)
cat(paste0("P-value for uniformity: ", pnorm(testUniformMST(bench=samples.rej, test=samples.sds)), "\n"))
