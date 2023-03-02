//clusterin.cpp
//implements clustering based on heatFlow and 1D clustering 

#include "clustering.h"


//selects start Node for heat dissipation to have high likelihood to lie in a cluster
Node* selectStartNode(const Graph* G,const std::unordered_map<Node*,int>* clusterLabels){
    //specifically selects Node that has the most neighbors
    //apply to kNN
    int maxNeighborCount=0;
    Node* startNode=nullptr;
    for (int i = 0; i < G->size(); i++)
    {
        Node* currNode=G->getVertex(i);
        if((*clusterLabels).at(currNode)!=-1)
        //this node is part of a previously detected cluster
        {
            continue;
        }
        int neighborCount=G->getVertex(i)->neighbors_.size();
        if(neighborCount>maxNeighborCount)
        {
            maxNeighborCount=neighborCount;
            startNode=G->getVertex(i);
        }    
    }
    return startNode;  
}


//performs heat clustering on data
//see clustering.h for detailed describtion of input/output format
std::vector<int> heatClustering(const std::vector<std::vector<double>>& data, double minClusterSize, double concentrationRadius){

    int k = std::log(data.size()+1)+1; //size of kNN to reduce computational cost in heat flow computations
    Graph* kNN = getMinkNN(data,k);
    std::unordered_map<Node*,int> clusterLabelsGraph = heatClustering(kNN, minClusterSize, concentrationRadius, true);
    int s=data.size();
    std::vector<int> clusterLabels(s,-1);
    for(const auto& [node, clusterLabel] : clusterLabelsGraph)
    {
        clusterLabels[kNN->getLabel(node)]=clusterLabel;
    }
    return clusterLabels;
}

//performs heat clustering on input data given as graph
//see clustering.h for detailed describtion of input/output format
std::unordered_map<Node*,int> heatClustering(Graph* G, double minClusterSize, double concentrationRadius, bool reduced){
    
    //GLOBAL PARAMETERS THAT MAY BE CHANGED HERE:

    //time used in one iteration of heat dissipation, smaller is more accurate but computationally more expensie
    double timeScale=0.01;
    //time we run heat dissipation before checking anything again, scaling by size or diam or.. of G is natural here
    double time=timeScale*G->size();
    //bandwidth parameter used in 1D clustering of heat vector via kernel density estimator
    double bandWidth=0.2*pow(G->size(),-0.2);
    //if maxDensity>=significance*minDensity in density estimate we separate into clusters
    double significance = 25; 
    //size of kNN to reduce computational cost in heat flow computations
    int k = std::log(G->size()+1)+1;
    //source is chosen below, if looking for a cluster containing a specific node it should be changed by hand
    
    //creat kNN graph
    Graph* kNN = new Graph(*G);
    if(!reduced){kNN = getMinkNN(G, k);}
   
    //replace euclidean distances by gaussian distances
    if(concentrationRadius!=0)
    {
        for (int i = 0; i < kNN->size(); i++)
        {
            Node* currNode = kNN->getVertex(i);
            auto it = currNode->neighbors_.begin();
            while(it!=currNode->neighbors_.end())
            {
                double dist = it->second;
                double expDist=std::pow(std::exp(-dist*concentrationRadius),-1);
                if(expDist==0)
                {
                    //erase edge
                    it=currNode->neighbors_.erase(it);
                }
                else
                {
                    currNode->neighbors_[it->first]=expDist/5; //done symmetrically anyways because we loop over all nodes
                    it++;
                }
            }
        }
    }

    //normalize distances in graph so that min distance is 1/2
    double minDist=std::numeric_limits<double>::infinity();    
    for(int i=0; i<kNN->size(); i++)
    {
        Node* currNode = kNN->getVertex(i);
        for(const auto& [nb, dist] : currNode->neighbors_)
        {
            minDist=std::min(minDist,dist);
        }
    }
    for(int i=0; i<kNN->size(); i++)
    {
        Node* currNode = kNN->getVertex(i);
        for(auto& [nb, dist] : currNode->neighbors_)
        {
            dist/=(2*minDist);
        }
    }
    
    std::unordered_map<Node*,int> clusterLabels;
    for(int i=0; i<G->size();i++)
    {
        clusterLabels[G->getVertex(i)]=-1;
    }
    
    std::unordered_map<Node*,int> clusterLabelskNN;
    for(int i=0; i<kNN->size();i++)
    {
        clusterLabelskNN[kNN->getVertex(i)]=-1;
    }

    //below is splitting into connected components and then running heatClustering on each of those
    if(!(kNN->isConnected()))
    {
        
       
        
        int numberOfComponents=kNN->connectedCompCount();
        std::vector<Graph*> components; //stores pointer to each component of kNN
        for (int i = 0; i < numberOfComponents; i++)
        {
            Graph* componentGraph = new Graph; //initialize empty
            components.push_back(componentGraph);
        }
        for (int i = 0; i < kNN->size(); i++)
        {
            Node* currNode = kNN->getVertex(i);
            int compLab = kNN->connectedCompLabel(currNode);
            components[compLab]->insert(currNode); //insert node into corresponding component
        }
        
        int clusterLabelCorrection=0; //on each component we label 0,..., need to correct for this globally
        for(int i=0; i<numberOfComponents; i++)
        {
            int temp=0; //stores max label occuring in this component
            if(components[i]->size()<minClusterSize*kNN->size()){continue;}
    
            double componentMinClusterSize=minClusterSize/(components[i]->size())*kNN->size();
            std::unordered_map<Node*,int> clusterLabelsOnComponent = heatClustering(components[i], componentMinClusterSize, 0); //defined on nodes of kNN
           
            for(const auto& [node, clusterLabel] : clusterLabelsOnComponent)
            {
                int idx = kNN->getLabel(node); //index of corresponding data point

                if(clusterLabel==-1)
                {
                    //case in which not assigned to any cluster within ith component
                    clusterLabels[G->getVertex(idx)]=-i-1;
                }
                else
                {
                    //case that it is assigned to a cluster within ith component
                    clusterLabels[G->getVertex(idx)]=clusterLabel+clusterLabelCorrection;
                }
                temp=std::max(temp,clusterLabel);
            }
            
            clusterLabelCorrection =  temp+1;
            std::cout << "Correction:" << clusterLabelCorrection << std::endl;
        }

        return clusterLabels;
    }
    
   std::cout << "size=" << kNN->size() <<std::endl;
  
    bool admissible=true; //have made sure ourselves above  
    double cutOffStep=minClusterSize/2; 
    int pointsLabeled=0; //keeps track of how many points we have identified to belong to clusters
    int count =0;
    double stop=1/(1-minClusterSize);

    int pointsLabeledCurrRound=-1;
    int startInd=0;
    int clusterLabelCorrection=0;



    while(pointsLabeled<(1-minClusterSize)*kNN->size())
    { 
       
        double cutOff=1/cutOffStep;
        //default choice is to select such that it is node with most neighbors,
        //high likelihood to sit in a cluster unless data highly symmetric
        //if initial source selection not successful try to use
        //node at idx=(minClusterSize/2*size)*i as sources
        Node* source;
        
        
        if(pointsLabeledCurrRound==0)
        {
            //last choice was not successfull
            //choose start node differently
            source = kNN->getVertex(startInd);
            startInd+=(minClusterSize*kNN->size()/2);
        }   
        else
        {
            source = selectStartNode(kNN, &clusterLabelskNN);
        }
        int sourceInd=kNN->getLabel(source);
        Node* sourceG = G->getVertex(sourceInd);

        //initialize heatDistribution to be concentrated in source
        //chosen so that equilibrium is 1 everyhwere
        std::unordered_map<Node*,double> heatDistribution;
        for(int i=0; i<kNN->size();i++)
        {
            Node* currNode=kNN->getVertex(i);
            if(currNode==source)
            {
                double heat=kNN->size();
                heatDistribution[currNode]=heat;
            }
            else
            {
                double cold=0;
                heatDistribution[currNode]=cold;
            }
        }

      

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
            
            if(heatLabel->size()>0)
            {
                //1D clustering ensures clusters has at least minCLusterSize
                int localCorrection=0;
                for(int i=0; i<kNN->size();i++)
                {
                    Node* currNodekNN=kNN->getVertex(i);
                    Node* currNodeG=G->getVertex(i);
                    auto it=(*heatLabel).find(heatDistribution[currNodekNN]);
                    if(it==(*heatLabel).end()){continue;}
                    else if(clusterLabels[currNodeG]==-1)
                    {
                        //checks if we have already labeled here
                        double heat = heatDistribution[currNodekNN];
                        if((*heatLabel)[heat]==-1){exit(123);}
                        clusterLabels[currNodeG]=(*heatLabel)[heat]+clusterLabelCorrection;
                        clusterLabelskNN[currNodekNN]=(*heatLabel)[heat]+clusterLabelCorrection;
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

    return clusterLabels;
}








