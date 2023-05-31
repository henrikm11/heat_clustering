//main.cpp

#include "clustering.h"
#include "graph.h"

template<typename Container>
void printVector(const Container& c){
    auto it=c.begin();
    auto stop=--c.end();
    std::cout<< "(";
    while(it!=stop){
        std::cout << *it << ",";
        it++;
    }
    std::cout << *it << ")" << std::endl;
    return;
}



int main(){

    std::vector<double> data = {0.01,0.02,0.03,0.04,0.002,0.001,0.056,0.2,-0.12,0.003,1.01,1.02,1.03,1.032,0.98,0.92,1.0123,0.4,0.41,0.43,0.425,0.38,0.35,0.39};
    std::vector<std::vector<double>> dataPoints(data.size());
    for(int i=0; i<data.size();i++){
        dataPoints[i]={data[i]};
    }
    

    std::vector<std::pair<int,double>> dataInd;
    for(size_t i=0; i<data.size(); i++){
        dataInd.push_back({i, data[i]});
    }

    double significance=2;
    double bandWidth = 0.02 ;
    double minClusterSize=0.05;

  
    std::unique_ptr<Graph> kNN = getkNN(dataPoints,7);
    if(kNN->checkSym()){
        std::cout <<"ok" <<std::endl;
    }
    else{
        std::cout << "not ok" << std::endl;
    }
    //bool connected=kNN->isConnected();
    //if(connected){std::cout << "connected" <<"\n";}
    //else{std::cout << "not connected" <<"\n";}



    //std::cout << kNN->vertices_[0]->coordinates_[0] << ",";
    
    
    for (size_t i = 0; i <kNN->size(); i++){
        std::cout << kNN->getVertex(i).getCoordinates()[0] << ",";
        //std::cout << kNN->getVertex(i).coordinates_[0] <<std::endl;
        //std::cout << count <<std::endl;
        //std::cout << "," << kNN->getVertex(i).neighbors_.size() << "\n";
    }
    
    //std::cout<<"here" <<std::endl;
   
    

    
    
    kNN->normalizeDistances(0.1);
    std::vector<double> initialDist(kNN->size(),0);
    initialDist[0]=1;

    
    
    printVector(initialDist);

   
    double timeScale=0.1;
    for(int i=0; i<300; i++){
        double high=-1;
        double low=-1;
        kNN->heatIterationStep(initialDist, timeScale, low, high);
       
    }

    printVector(initialDist);

    std::vector<int> clusterLabels=oneDimensionalClustering(initialDist, significance, bandWidth, minClusterSize);
    printVector(clusterLabels);
    //std::vector<int> clusterLabels=oneDimensionalClustering(data, significance, bandWidth, minClusterSize);
   
    //printVector(clusterLabels);

    return 0;



    //std::cout << baseDensity << std::endl;
    //std::cout << expectedCount << std::endl;

    /*
    OneDimClusterHelper Helper = OneDimClusterHelper(); //private again
    double baseDensity = Helper.getBaseDensity(dataInd);
    double expectedCount = Helper.getExpectedBaseWindowCount(dataInd, bandWidth);
    std::vector<std::pair<int,int>> clusterSplits = Helper.getClusterSplits(dataInd, significance, bandWidth, minClusterSize);
    

    //std::cout << densities.first << "\n" << densities.second << "\n";
    //std::cout << clusterSplits.size();

    for(const auto& split : clusterSplits){
        std::cout << split.first << "," << split.second << std::endl;
    }

    return 0;
    


    std::vector<int> labels = oneDimensionalClustering(data, 1.1, 0.1, 0.1);
    //data, significance, bandWidth, minClusterSize
    printIntVector(labels);
   
    return 0;
    */


    



    /*
    Node* origin = new Node(std::vector<double>(1,0)); //node at the 1D origin
    Node* nextToOrigin = new Node(std::vector<double>(1,1));
    Node* farAway = new Node(std::vector<double>(1,100));
    //Node* veryfarAway = new Node(std::vector<double>(1,1000));
    origin->neighbors_.push_back({nextToOrigin,1});
    origin->neighbors_.push_back({farAway,100});
    nextToOrigin->neighbors_.push_back({farAway,99});
    //farAway->neighbors_.push_back({veryfarAway,999});

    std::vector<Node*> vertices={origin, nextToOrigin, farAway};

    Graph G = Graph(vertices);

    std::cout<< G.getVertex(0).neighbors_.size();
    G.remove(&G.getVertex(2));
    std::cout<< G.getVertex(0).neighbors_.size();
    */



    //std::cout << G.getVertex(2).neighbors_.size();
    
    /*
    std::vector<double> initialDist = {0,1,0};
    double timeScale=1e-3;
    for (int i = 0; i < 10000; i++)
    {
        G.heatIterationStep(initialDist, timeScale);
        if(i%100==0){
            std::cout<<initialDist[0] << "," << initialDist[1] << "," << initialDist[2];
            std::cout << std::endl;
        }

    }
    */

}