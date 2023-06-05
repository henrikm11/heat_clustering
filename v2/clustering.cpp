//clustering.cpp

//implementation of heat clustering algorithm

#include "graph.h"
#include "clustering.h"

//helper class to split heatClustering function into various methods not to be accessed elsewhere
class ClusteringHelper{
private:
    friend std::vector<int> heatClustering(const Graph& G, double minClusterSize, double concentrationRadius, double significance, double bandWidth, double timeScale, bool reduced);
    int selectStartNode(const Graph& G, const std::vector<int>& clusterLabels);
    std::vector<int> heatClusteringConnected(const Graph& G, double minClusterSize, double concentrationRadius, double significance, double bandWidth, double timeScale);
};

int ClusteringHelper::selectStartNode(const Graph& G, const std::vector<int>& clusterLabels){
    //specifically selects Node that has the most neighbors
    //apply to kNN
    int maxNeighborCount=0;
    int startNodeIdx=-1;
    for (int i = 0; i < G.size(); i++)
    {
        if(clusterLabels[i]!=-1){continue;}
        //this node is part of a previously detected cluster
      
        int neighborCount=G.getVertex(i).neighborCount();
        if(neighborCount>maxNeighborCount)
        {
            maxNeighborCount=neighborCount;
            startNodeIdx=i;
        }    
    }
    return startNodeIdx;  
}

std::vector<int> heatClustering(const std::vector<std::vector<double>>& data,double minClusterSize, double concentrationRadius, double significance, double bandWidth, double timeScale){
    int k = std::log(data.size()+1)+1;
    std::unique_ptr<Graph> kNN=getkNN(data,k);
    return heatClustering(*kNN, minClusterSize, concentrationRadius, significance, bandWidth, timeScale, true);
};

std::vector<int> heatClustering(const Graph& G, double minClusterSize, double concentrationRadius, double significance, double bandWidth, double timeScale, bool reduced=false){
    std::vector<int> clusterLabels(G.size(),-1);
    int k = std::log(G.size()+1)+1;
    std::unique_ptr<Graph> kNN=getkNN(G,k);
    ClusteringHelper Helper;
    std::vector<int> componentLabels=kNN->getComponentLabels();

    /*
    
    TO DO
    
    */

   //split into components
   //apply heatClusteringConnected to components
   //collect labels from those and update clusterLabesl

    return clusterLabels;
}

std::vector<int> ClusteringHelper::heatClusteringConnected(const Graph& G, double minClusterSize, double concentrationRadius, double significance, double bandWidth, double timeScale){
    std::vector<int> clusterLabels(G.size(),-1);

    return clusterLabels;
}



