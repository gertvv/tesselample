source('fast_sampling.R')

## n-dimensional
for (n in 3:9)  {
constr <- do.call(mergeConstraints,
  c(list(simplexConstraints(n=n)),
    lapply(2:n, function(i) { ordinalConstraint(n=n, i=i, j=1) })))

  transform <- simplex.createTransform(n=n)
  constr <- transformConstraints(transform, filterConstraints(constr, constr$dir == "<="))

  result <- sds(constr, 1E4, homogeneous=TRUE)
  print(paste("n =", n))
  print(result$timing)
}
