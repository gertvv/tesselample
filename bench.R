source('fast_sampling.R')

## n-dimensional
for (n in 3:10)  {
constr <- do.call(mergeConstraints,
  c(list(simplexConstraints(n=n)),
    lapply(2:n, function(i) { ordinalConstraint(n=n, i=i, j=1) })))

  transform <- simplex.createTransform(n=n)
  constr <- transformConstraints(transform, filterConstraints(constr, constr$dir == "<="))

  time <- system.time(
    result <- sds(constr, 1E4, homogeneous=TRUE)
  )
  print(paste("n + 1 =", n))
  print(result$timing)
  print("=== Total ===")
  print(time)
}
