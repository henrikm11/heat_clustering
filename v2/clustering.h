//clustering.h
//single header file for clustering project
//included in
//main.cpp
//graph.cpp

#ifndef CLUSTERING_H
#define CLUSTERING_H

#include "graph.h"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <exception>
#include <memory>
#include <queue>
#include <algorithm>


/// @brief detects clusters in one dimensional data via kernel density estimate
/// @param data input data
/// @param significance quotient of densities to mark split between clusters
/// @param bandWidth size of windows used to estimate density
/// @param minClusterSize relative size minimally required to be a cluster
/// @return vector of labels, indexed as input data
std::vector<int> oneDimensionalClustering(const std::vector<double>& data, const double significance, const double bandWidth, const double minClusterSize);




std::vector<int> heatClustering(const Graph& G, double minClusterSize, double concentrationRadius, double significance, double bandWidth, double timeScale, bool reduced=false);
std::vector<int> heatClustering(const std::vector<std::vector<double>>& data,double minClusterSize, double concentrationRadius, double significance, double bandWidth, double timeScale);


//Helper class to implement oneDimensionalClustering
class OneDimClusterHelper{
public: //change back to private
    friend std::vector<int> oneDimensionalClustering(const std::vector<double>&, const double, const double, const double);

    /// @brief returns density in bandWidth windows if data is uniformly distributed up to some noise
    /// @param dataInd data with indices, sorted by data
    /// @return -1 if dataInd.size==0 or all data at a single point up to noise
    double getBaseDensity(const std::vector<std::pair<int,double>>& dataInd);

    /// @brief returns expected number of data points in window of size bandWidth if data is uniformly distributed up to some noise
    /// @param dataInd data with indices, sorted by data
    /// @param bandWidth size of sampling window
    /// @return expected count or dataInd.size() if points essentially concentrated at one point
    double getExpectedBaseWindowCount(const std::vector<std::pair<int,double>>& dataInd, const double bandWidth);

    std::vector<std::pair<int,int>> getClusterSplits(const std::vector<std::pair<int,double>>& dataInd, const double significance, const double bandWidth, const double minClusterSize);
};

#endif //CLUSTERING_H
