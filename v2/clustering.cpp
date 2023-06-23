//clustering.cpp

//implementation of heat clustering algorithm

#include "graph.h"
#include "clustering.h"




std::vector<int> heatClustering(
    const std::vector<std::vector<double>>& data,double minClusterSize,
    double concentrationRadius,
    double significance, 
    double bandWidth, 
    double timeScale
    )
{
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
    )
{
    std::vector<int> clusterLabels(G.size(),-1); //return variable

    int k = std::log(G.size()+1)+1; //size of kNN 
    std::unique_ptr<Graph> kNN=reduce(G,k,reduced); //get kNN
    
    //split kNN into components, smart pointer does memory management for us
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

    //get labels on components using Helper function
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

    //correction to create global labels -  labels on components are always 0,1,2,...
    std::vector<int> labelCorrections(componentCount,0);
    int correction=0;
    for(int i=0; i<componentCount; i++){
        labelCorrections[i]=correction;
        correction+=*(std::max(clusterLabelsOnComponents[i].begin(), clusterLabelsOnComponents[i].end()));
        correction++;
    }

    //create global labels from labels on components
    std::vector<int> posInComponent(componentCount,0);
    for(size_t i=0; i<kNN->size(); i++){
        int currComponent=componentLabels[i];
        if(clusterLabelsOnComponents[currComponent][posInComponent[currComponent]]!=-1){
            //otherwise not identified as part of cluster
            clusterLabels[i]=clusterLabelsOnComponents[currComponent][posInComponent[currComponent]]+labelCorrections[currComponent];
        }
        posInComponent[currComponent]++;
    }

    return clusterLabels;
}


int ClusteringHelper::selectSourceNode(
    const Graph& G,
    const int minClusterSize,
    const std::vector<int>& clusterLabels,
    int& startInd,
    const int pointsLabeledCurrRound
    )
{
    //selects Node that has the most neighbors
    //if heat flow from that node does not help we select a somewhat random nodoe to try instea
    //returns -1 if no good choices left
    int sourceNodeIdx=-1; //return value

    if(pointsLabeledCurrRound==0)
        //choice of source wasn't good, choose random node
    {
        //last choice was not successfull
        //choose start node differently
        sourceNodeIdx=startInd;
        startInd+=(minClusterSize*G.size()/2);
        if(sourceNodeIdx>=G.size()){sourceNodeIdx=-1;} //check if all hope is lost
        return sourceNodeIdx;  
    }  
    
    
    int maxNeighborCount=0;
    for (int i = 0; i < G.size(); i++)
    {
        if(clusterLabels[i]!=-1){continue;}
        //this node is part of a previously detected cluster
    
        int neighborCount=G.getVertex(i).neighborCount();
        if(neighborCount>maxNeighborCount)
        {
            maxNeighborCount=neighborCount;
            sourceNodeIdx=i;
        }    
    }
    return sourceNodeIdx;  
}

std::vector<int> ClusteringHelper::heatClusteringConnected(
    Graph& G, 
    double minClusterSize, 
    double concentrationRadius, 
    double significance, 
    double bandWidth, 
    double timeScale
    )
{
    //knows that G is kNN
    assert(G.isConnected());
    G.gaussianDistances(concentrationRadius); //change distances to be Gaussian and normalizes
    std::vector<int> clusterLabels(G.size(),-1);


    double cutOffStep=minClusterSize/2; 
    int pointsLabeled=0; //keeps track of how many points we have identified to belong to clusters
    int count = 0;
    double stop=1/(1-minClusterSize);

    int pointsLabeledCurrRound=-1; //initialized to -1 to not skip first round of smart source selection
    int startInd=0;
    int clusterLabelCorrection=0; //OneDim clustering labels always start at 0, may need to correct for that


    while(pointsLabeled<(1-minClusterSize)*G.size()){
        double cutOff=1/cutOffStep;


        //select node to heat up initially
        int sourceNodeIdx=selectSourceNode(G, minClusterSize, clusterLabels,startInd,pointsLabeledCurrRound);
        if(sourceNodeIdx==-1){break;} //all hope to find clusters is lost

        std::vector<double> heatDistribution(G.size(),0);
        heatDistribution[sourceNodeIdx]=1;
        //equilibrium is 1/G.size() everywhere

        pointsLabeledCurrRound=0;
        //running heat from source and attempting clustering
        while(cutOff>stop){
        //if cluster containing source contains x*size of all points 
        //we expect it to find when heatDistribution[source]~1/(x*size)
            const double inf = std::numeric_limits<double>::infinity();
            double heatSourcePrev = inf, gradientHeatSource=inf, maxGradient=-1;
            double low=0, high=1; //min, max of current heat Distribution, updated by heatIteration
            
            while(heatDistribution[sourceNodeIdx]>cutOff){//let heat dissipate 
                G.heatIterationStep(heatDistribution, timeScale, low, high);
                gradientHeatSource = std::abs(heatDistribution[sourceNodeIdx]-heatSourcePrev);
                if(gradientHeatSource<inf){maxGradient = std::max(maxGradient, gradientHeatSource);}
                if((significance)*(significance)*gradientHeatSource < maxGradient){
                    //small (time) gradient at source, attempt 1D clustering
                    cutOffStep-=minClusterSize; //this way the cuttOff is changed back to same value below
                    break;
                }          
                heatSourcePrev = heatDistribution[sourceNodeIdx];
            }

            //we can now use 1D clustering to check if we have already found clusters
            // otherwise we keep dissipating heat

            std::vector<int> heatLabels=oneDimensionalClustering(
                heatDistribution,
                significance,
                bandWidth,
                minClusterSize
            );
            
            bool foundCluster=false;
            for(int label : heatLabels){ //check if cluster found
                if(label!=-1){
                    foundCluster=true;
                    break;
                }
            }
            if(foundCluster){ //label points accordingly
                //1D clustering ensures clusters has at least minCLusterSize
                int localCorrection=0;
                for(int i=0; i<G.size();i++){
                    if(heatLabels[i]!=-1 && clusterLabels[i]==-1){
                        clusterLabels[i]=heatLabels[i]+clusterLabelCorrection;
                        pointsLabeledCurrRound++;
                        pointsLabeled++;
                        localCorrection=std::max(localCorrection, heatLabels[i]);
                    }
                }
                clusterLabelCorrection+=localCorrection+1;
                //more efficient so restart heatflow from another source, hence we stop here
                break; 
            }
            
    /*
    TO DO
    finish adjustment of old code
    */




        }


    }



    /*
    { 
       
    

        pointsLabeledCurrRound=0;
        //running heat from source and attempting clustering
        while(cutOff>stop)
        {
        //if cluster containing source contains x*size of all points we expect it to find when heatDistribution[source]~1/x
            const double inf = std::numeric_limits<double>::infinity();
            double heatSourcePrev = inf;
            double gradientHeatSource;
            double maxGradient = -1;
        
            while(heatDistribution[source]>cutOff)
            {
                //kNN->heatIterationStep(timeScale, heatDistribution, admissible);
                kNN->heatDiffusion(time, timeScale, heatDistribution,admissible);
                gradientHeatSource = std::abs(heatDistribution[source]-heatSourcePrev);
                if(gradientHeatSource<inf){maxGradient = std::max(maxGradient, gradientHeatSource);}
                if((significance)*(significance)*gradientHeatSource < maxGradient)
                {
                    cutOffStep-=minClusterSize; //this way the cuttOff is changed back to same value below
                    break;
                }
                
                heatSourcePrev = heatDistribution[source];
                count++;
            }
           
        
            //we can now use 1D clustering to check if we have already found clusters
            // otherwise we keep dissipating heat
            std::vector<double> heatVector;
            for(int i=0; i<kNN->size(); i++)
            {
                Node* currNode = kNN->getVertex(i);
                heatVector.push_back(heatDistribution[currNode]);
            }
        
            //1D clustering
            std::unordered_map<double,int>* heatLabel=oneDimensionalClustering(heatVector, significance, bandWidth, minClusterSize);
            

            //HERE


            if(heatLabel->size()>0)
            {
                //1D clustering ensures clusters has at least minCLusterSize
                int localCorrection=0;
                for(int i=0; i<kNN->size();i++)
                {
                    Node* currNode=kNN->getVertex(i);
                    auto it=(*heatLabel).find(heatDistribution[currNode]);
                    if(it==(*heatLabel).end()){continue;}
                    else if(clusterLabels[currNode]==-1)
                    {
                        //checks if we have already labeled here
                        double heat = heatDistribution[currNode];
                        if((*heatLabel)[heat]==-1){exit(123);}
                        clusterLabels[currNode]=(*heatLabel)[heat]+clusterLabelCorrection;
                        pointsLabeledCurrRound++;
                        pointsLabeled++;
                        localCorrection=std::max(localCorrection,(*heatLabel)[heat]);
                    }
                }
                clusterLabelCorrection+=localCorrection+1;
                //more efficient so restart heatflow from another source, hence we stop here
                break; 
            }
            cutOffStep+=minClusterSize;
            cutOff=1/cutOffStep;
            if(pointsLabeled>(1-minClusterSize)*kNN->size()){break;}
        }
        if(startInd>=kNN->size()){break;}
    }   

    */


    return clusterLabels;
}



