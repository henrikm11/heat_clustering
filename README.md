# Heat-clustering

C++ code for clustering using heat flow

## Some background
The algorithm first constructs a graph from a vector of data using a version of k nearest neighbors. Then heat dissipation on this graph is used as a dimension reduction map. The resulting one dimensional data points are clustered using a kernel density estimator.

## Parameters
-) ConcentrationRadius, see kNN Graph
-) timeScale, see heat dissipation
-) significance, see 1D clustering

## Implementation

Some details on key elements of the implementation

### k nearest neighbor graph
Construct symmetric kNN graph, default value for k is roughly k~log(n), where n is size of data set. Distance are computed either with respect to Euclidean distance or using Gaussian kernel (normalized to have max 1) with mean 0 and variance conentrationRadius/2.

### heat dissipation
Heat dissipation is modeled using a simple Euler approximation for the heat equation, timeScale is the size of time steps in approximation.
Initial conditions are 0 everywhere except in a source Node.
Source node is selected to have high likelihood to be contained in a cluster.

### dimension reduction
Once a source and a time are selected the Graph is mapped to a one dimensional vector by checking heat at nodes at that time given intial distribution described above.
We do so whenever gradient (in time) of heat becomes small at source node.

### 1D clustering
This is using step functions of fixed bandwidth depending on data size to approximate the density function. We label a cluster if density between points of high density drops significantly. This is tuned using parameter significance. Higher significance means clusters need to be separated more clearly.
