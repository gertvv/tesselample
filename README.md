# Uniform Sampling from Tesselations of a Convex Polytope

This is some simple R code that samples uniformly from convex polytopes. It does it in the following way:

 - Find the vertices of the polytope (its convex hull)
 - Tesselate (triangulate) the set of vertices, resulting in simplices
 - Draw a uniform value to decide which simplex to sample from
 - Sample a point uniformly from the simplex

This works reasonably well up to 9 dimensions or so, at which point the number of simplices is [simply too large](http://math.stackexchange.com/questions/474857/find-the-smallest-triangulation-of-the-n-dimensional), growing as n!, making the tesselation impossible.

Further documentation, including proofs and timing results, is available in the LaTeX document `proof.tex`. The implementation is in `fast_sampling.R` and a simple test script in `test_sampling.R`. A simple benchmark is in `bench.R`.
