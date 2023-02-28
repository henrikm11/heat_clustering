# Heat-clustering

C++ code for clustering using heat flow

## Some background
The algorithm first constructs a graph from a vector of data using a version of k nearest neighbors. Then heat dissipation on this graph is used as a dimension reduction map. The resulting one dimensional data points are clustered using a kernel density estimator.

## Implementation Details

Details about the implementation of clustering algorithm
