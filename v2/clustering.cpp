//clustering.cpp

//implementation of heat clustering algorithm

#include "graph.h"
#include "clustering.h"




std::vector<int> heatClustering(const std::vector<std::vector<double>>& data,double minClusterSize, double concentrationRadius, double significance, double bandWidth, double timeScale){
    int k = std::log(data.size()+1)+1;
    std::unique_ptr<Graph> kNN=getkNN(data,k);
    return heatClustering(*kNN, minClusterSize, concentrationRadius, significance, bandWidth, timeScale, true);
};

std::unique_ptr<Graph> reduce(const Graph& G, int k, bool reduced){
    if(reduced){
        //already reduced, simply copy
        std::unique_ptr<Graph> reducedGraph(new Graph(G));
        return reducedGraph;
    }
    std::unique_ptr<Graph> kNN = getkNN(G,k);
    return kNN;
}

std::vector<int> heatClustering(
    const Graph& G, 
    double minClusterSize, 
    double concentrationRadius, 
    double significance, 
    double bandWidth, 
    double timeScale, 
    bool reduced
    ){
    std::vector<int> clusterLabels(G.size(),-1); //return variable
    int k = std::log(G.size()+1)+1; //size of kNN 
    std::unique_ptr<Graph> kNN=reduce(G,k,reduced); //get kNN
    
    //split kNN into components
    std::vector<int> componentLabels=kNN->getComponentLabels();
    int componentCount =  kNN->componentCount();
    std::vector<std::vector<Node*>> componentVertices(componentCount);
    for(int i=0; i<kNN->size(); i++){
        componentVertices[componentLabels[i]].push_back(&(kNN->getVertex(i)));
    }
    std::vector<std::unique_ptr<Graph>> components(componentCount);
    for(int j=0; j<componentCount; j++){
        std::unique_ptr<Graph> component(new Graph(componentVertices[j]));
        components[j]=std::move(component);
    }

    //get labels on components
    ClusteringHelper Helper;
    std::vector<std::vector<int>> clusterLabelsOnComponents;
    for(size_t i=0; i<componentCount; i++){
        double correctedMinClusterSize=minClusterSize*G.size()/components[i]->size(); //clusters of this size in components[i] are significant in all of G
        
        std::vector<int> clusterLabelsComponent = Helper.heatClusteringConnected(*components[i],
            correctedMinClusterSize, 
            concentrationRadius, 
            significance, 
            bandWidth,
            timeScale
        );
        clusterLabelsOnComponents.push_back(clusterLabelsComponent);
    }


    /*
    
    TO DO
    
    */
   //update clusterLabes globally using labels on components

    return clusterLabels;
}


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

std::vector<int> ClusteringHelper::heatClusteringConnected(const Graph& G, double minClusterSize, double concentrationRadius, double significance, double bandWidth, double timeScale){
    assert(G.isConnected());
    std::vector<int> clusterLabels(G.size(),-1);


    return clusterLabels;
}



