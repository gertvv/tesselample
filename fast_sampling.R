## TA Trikalinos and GHM van Valkenhoef, Nov 2014 
## Efficient sampling from uniform densities over n-polytopes
library(geometry)
library(hitandrun)

## A d-polytope P can be decomposed into d-dimensional simplices S_1, ..., S_K,
## such that the union of S_k (k=1,...K) equals P, and their intersection is
## the empty set or a (d-1)-dimensional simplex. 
## A Delaunay triangulation is such a simplicial decomposition that maximizes
## the minimum angle of each S_k. 
## The problem then reduces to uniform sampling from these simplices.
sds <- function(constr, N, homogeneous=FALSE, transform=NULL) {
  timing <- list()
 
  if (homogeneous) { # eliminate the homogeneous coordinate
    n <- ncol(constr$constr)
    constr$rhs <- constr$rhs - constr$constr[, n]
    constr$constr <- constr$constr[,1:(n-1)]
  }
  timing$vertices <- system.time(
    verts <- findVertices(constr)
  )
  n <- ncol(verts)
  timing$tesselation <- system.time(
    decomposition <- delaunayn(verts, options="QJ", full=TRUE)
  )
  simplex.vert.indices <- decomposition$tri

  ## the number of simplices is K
  K <- nrow(simplex.vert.indices)
  volumes <- decomposition$areas
  prob <- volumes / sum(volumes)

  timing$info <- c(n=n, v=nrow(verts), K=K)

  ## for each simplex k, do a convex combination of it's vertices
  timing$sample.indices <- system.time(
    simplex.idx <- sample(K, N, replace=TRUE, prob=prob)
  )
  samples <- matrix(NA, nrow=N, ncol=n)
  timing$sample.weights <- system.time(
    W <- simplex.sample(n + 1, N)$samples
  )

  timing$transform <- system.time(
    for(i in 1:N){
      # the k-th simplex 
      S.k.prime <- verts[simplex.vert.indices[simplex.idx[i],],]
      
      # random uniform convex combination weights
      samples[i, ] <- W[i, , drop=FALSE] %*% S.k.prime
    }
  )

  # Add the homogeneous coordinate if necessary
  samples <- 
    if (homogeneous) {
      cbind(samples, 1)
    } else {
      samples
    }

  # Apply the transformation if specified
  if (is.null(transform)) {
    list(samples=samples,timing=timing)
  } else {
    list(samples=samples %*% t(transform), timing=timing)
  }
}

# filter a set of constraints
filterConstraints <- function(constr, sel) {
  list(constr = constr[['constr']][sel, , drop=FALSE],
       rhs = constr[['rhs']][sel],
       dir = constr[['dir']][sel])
}
