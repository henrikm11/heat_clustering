//knn.cpp

#include "clustering.h"
#include "graph.h"


//only use in implementation of getkNN
//leakage easily occurs if used by itself
class kNNHelper{
private:
    Graph* getPrekNN(std::vector<std::vector<double>>& data, int k);
    Graph* getPrekNN(const Graph&, int k); //creates new graph
    Graph* reducePrekNN(Graph* prekNN); //removes directed edges wo/ counterpart from preKNN
    std::vector<std::vector<std::pair<int,double>>> computeDistances(const std::vector<std::vector<double>>& data);
    std::vector<int> getMults(std::vector<std::vector<std::pair<int,double>>>& distances);
    void computeAdjLists(Graph& G);
    void removeNullity(Graph& prekNN);
    friend std::unique_ptr<Graph> getkNN(std::vector<std::vector<double>>& data, int k);
    friend std::unique_ptr<Graph> getkNN(const Graph& G, int k);
};


std::vector<std::vector<std::pair<int,double>>> kNNHelper::computeDistances(const std::vector<std::vector<double>>& data){
    //distances[i][j]={j,dist[i][j]} for later use
    std::pair<int,double> initPair(0,0);
    std::vector<std::pair<int,double>> initVec(data.size(), initPair);
    std::vector<std::vector<std::pair<int,double>>> distances(data.size(),initVec);
    
    if(data.size()==0){return distances;}
    int d = data[0].size();
    for (size_t i = 0; i < data.size(); i++){  
        for (size_t j = 0; j<=i; j++){  
            //j=i case required to correct index in pair!
            double dist=0; 
            for(size_t l = 0; l < d; l++){   
                //square of standard Euclidean distances
                dist+=((data[i][l]-data[j][l])*(data[i][l]-data[j][l]));
            }
            distances[i][j]=std::pair<int,double>({j,dist});
            distances[j][i]=std::pair<int,double>({i,dist});
        }
    }
    return distances;
}

/// @brief computes multiplicity of data points, only works if ordering is the same in each column
/// @param distances matrix of distances
/// @return vector giving multiplicities shifted to first occurence
std::vector<int> kNNHelper::getMults(std::vector<std::vector<std::pair<int,double>>>& distances){
    std::vector<int> mults(distances.size(),1);
    for(int i=distances.size()-1; i>-1; i--){
        for(int j=distances[i].size()-1; j>i;j--){
            auto nb=distances[i][j]; //nb.first=j
            if(nb.second!=0){continue;}
            mults[i]+=mults[j];
            mults[j]=0;
        }
    }
    return mults;
}

void kNNHelper::computeAdjLists(Graph& G){
    std::unordered_map<Node*,int> labels;
    for(size_t i=0; i<G.size(); i++){
        labels[G.vertices_[i]]=i;
    }
    G.adjacencyLists_.clear();
    for(size_t i=0; i<G.size(); i++){
        std::vector<std::pair<int,double>> adj;
        for(const auto& edge : G.vertices_[i]->neighbors_){
            adj.push_back({labels[edge.first],edge.second});
        }
        G.adjacencyLists_.push_back(adj);
    }
    return; 
}



/// @brief returns directed kNN from data, designed for case k<=data.size()/2
/// @param data vector of input data
/// @param k number of neighbors
/// @return directed graph of k nearest neighbors
Graph* kNNHelper::getPrekNN(std::vector<std::vector<double>>& data, int k){
    assert(k<=data.size()-1);
    Graph* prekNN = new Graph();
    if(data.size()==0){return prekNN;} //no data
   
    //compute distances and multiplicities 
    std::vector<std::vector<std::pair<int,double>>> distances = computeDistances(data);
    std::vector<int> mults=getMults(distances); //mults[i]=multiplicity of node i, redundancies removed

    //add nodes with correct multiplicity to prekNN
    prekNN->vertices_.reserve(data.size());
    std::vector<int> corrections(data.size(),0); //data[i]<->vertices[i-corrections[i]] if mults[i]>0
    int correction=0;
    for(size_t i=0; i<mults.size();i++){
        if(mults[i]!=0){prekNN->vertices_.push_back(new Node(data[i],mults[i]));}
        else{correction++;}
        corrections[i]=correction; 
    }



    //add neighbors
    auto cmp = [] (const std::pair<int,double>& a, const std::pair<int,double>& b){return a.second<b.second;};
    std::vector<int> counts(prekNN->size(),0);
    for(int i=0; i<data.size(); i++){
        //find k nearest neighbors of data[i]
        if(mults[i]==0){continue;} //data point is a duplicate
        
        int idx1 = i-corrections[i]; //vertices_[idx1] is the corresponding node
        std::nth_element(distances[idx1].begin(), distances[idx1].begin()+k, distances[idx1].end(),cmp);
        double threshold = distances[idx1][k].second; //check if any more elements of distance threshold in list
        std::vector<int> counts(prekNN->size(),0);
        for(const auto& nb : distances[idx1]){
            int j = nb.first;
            if(mults[j]==0){continue;} //edges to non existent nodes
            int idx2 = j-corrections[j];
            if(idx1==idx2){continue;} //no self edges
            double dist = nb.second;
            if(dist<=threshold){
                std::pair<Node*,double> edge={prekNN->vertices_[idx2],dist};
                prekNN->vertices_[idx1]->addNeighbor(edge);
            }
        }
    }
   

    return prekNN;
}


void kNNHelper::removeNullity(Graph& prekNN){
    size_t size=prekNN.size();
    size_t i=0;
    while(i<size){
        if(prekNN.vertices_[i]->mult_==0){
            prekNN.remove(prekNN.vertices_[i]);
            size--;
        }
        else{i++;}
    }
    return;
}

/// @brief TODO!!!
Graph* kNNHelper::getPrekNN(const Graph& G, int k){
    Graph* prekNN = new Graph(G);
    //removeNullity(*prekNN);
    
    //sort adjacencyLists and remove neighbors
    auto cmp = [](const std::pair<int,double> a, const std::pair<int,double> b){return a.second<b.second;};
    for(size_t i=0; i<prekNN->size(); i++){
        if(prekNN->adjacencyLists_[i].size()<k){continue;} //nothing to be removes
        std::nth_element(prekNN->adjacencyLists_[i].begin(),
            prekNN->adjacencyLists_[i].begin()+k-1,
            prekNN->adjacencyLists_[i].end(),
            cmp
        );
        double threshold =prekNN->adjacencyLists_[i][k-1].second;

        //now delete everything above threshold
        for(auto it = prekNN->adjacencyLists_[i].begin();it!=prekNN->adjacencyLists_[i].begin();it++){
            if(it->second<=threshold){continue;}
            prekNN->vertices_[i]->removeNeighbor(*(prekNN->vertices_[it->first]));
        }
    }
    computeAdjLists(*prekNN);
    return prekNN;
}


template<class T>
bool checkSymMat(std::vector<std::vector<T>> mat){
    int n=mat.size();
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(mat[i][j]!=mat[j][i]){return false;}
        }
    }
    return true;
}


Graph* kNNHelper::reducePrekNN(Graph* G){
    //reduce prekNN to kNN  
    //pass on ownership
    if(G==nullptr){return G;}
    std::vector<std::vector<int>> adjMat(G->size(),std::vector<int>(G->size(),0)); //adjMat[i][j]=1 iff there is an edge from i to j
    //adjMat==0 checked
    std::unordered_map<Node*,int> labels;
    for(int i=0; i<G->size(); i++){
        labels[G->vertices_[i]]=i;
    }
    for(int i=0; i<G->size(); i++){
        for(const auto& edge : G->vertices_[i]->neighbors_){
            adjMat[i][labels[edge.first]]=1; //edge from i to labels[...]
        }
    }

    for(int i=0; i<G->size(); i++){
        std::vector<int> counts(G->size(),0);
        for(auto& edge : G->vertices_[i]->neighbors_){
            int j = labels[edge.first];
            counts[j]++;
            if(counts[j]>1){std::cout <<"confused";}
        }
    }

    for(int i=0; i<G->size(); i++){
        //std::vector<int> counts(G->size(),0);
        for(auto it=G->vertices_[i]->neighbors_.begin(); it!=G->vertices_[i]->neighbors_.end();){
            std::pair<Node*,double> edge=*it;
            int j = labels[edge.first];
            //counts[j]++;
            //if(counts[j]>1){std::cout <<"confused here";}
            if(adjMat[j][i]==0){
                it=G->vertices_[i]->removeNeighborIt(*(edge.first));
                adjMat[i][j]=0;
            }
            else{
                it++;
            }
        }
    }
    assert(checkSymMat(adjMat));
    return G;    
}

std::unique_ptr<Graph> getkNN(std::vector<std::vector<double>>& data, int k){
    kNNHelper Helper;
    Graph* kNN = Helper.reducePrekNN(Helper.getPrekNN(data,k));
    Helper.computeAdjLists(*kNN);
    std::unique_ptr<Graph> kNNPtr(kNN); //prevents leakage
    return kNNPtr;
}


std::unique_ptr<Graph> getkNN(const Graph& G, int k){
    kNNHelper Helper;
    Graph* kNN = Helper.reducePrekNN(Helper.getPrekNN(G,k));
    Helper.computeAdjLists(*kNN);
    std::unique_ptr<Graph> kNNPtr(kNN); //prevents leakage
    return kNNPtr;
}
