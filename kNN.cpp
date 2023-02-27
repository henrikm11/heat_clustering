//kNN.cpp

//computes various types of kNN graphs for data set
//algorithm is basically brute force which is non optimal in low dimensions
//only optimization is use of priority queues

//in all functions the order of nodes/labeling is preserved!

//as of now it only uses square of Euclidean distance as weights in graph

#include "clustering.h"


//returns directed kNN, i.e. there is an edge from a to b iff b is among k closest to a
Graph* getPrekNN(const std::vector<std::vector<double>>& data, int k){
    //only used to efficiently computed actual kNNs later
    if(data.size()==0){
        Graph* kNN = new Graph;
        return kNN;     
    }

    //compute distances and store them in nxn matrix distances
    int n=data.size();
    int d=data[0].size(); //dimension of input data
    std::vector<double>temp(n,0);
    std::vector<std::vector<double>> distances(n,temp); //initialized to all 0

    for (int i = 0; i < n; i++)
    {  
        for (int j = 0; j < i; j++)
        {   
            for(int l = 0; l < d; l++)
            {   
                //so far we use (square of) standard Euclidean distances
                distances[i][j]+=((data[i][l]-data[j][l])*(data[i][l]-data[j][l]));
            }
           
            distances[j][i]=distances[i][j];
        }
    }

    //create Nodes of graph with no neighbors yet
    std::vector<Node*> Nodes;
    for (int i=0; i<data.size(); i++)
    {
        Node* temp = new Node(data[i]);
        Nodes.push_back(temp);
    }
    //coordinates associated via Nodes[i]<->data[i]

    //handle case of k>=n-1
    if(data.size()<=k-1){
        //put neighbors to all the nodes
        Graph* kNN = new Graph(Nodes);
        for(int j=0; j<data.size(); j++){
            for(int idx=0; idx<data.size(); idx++){
                if(idx==j){continue;}
                Nodes[j]->neighbors_[Nodes[idx]]=distances[j][idx]; //only add this direction, not symmetrically!
            }
        }
        return kNN;
    }
      
    //add k nearest neighbors to each Node
    for(int i=0; i<data.size(); i++){
        //compute k nearest neighbors of data[i]

        auto cmp = [](std::pair<int,double> a, std::pair<int,double> b){
            return (a.second<b.second); //this makes Q max heap
        };
        std::priority_queue<std::pair<int,double>, std::vector<std::pair<int,double>>, decltype(cmp)> Q(cmp);

        bool duplicate=false; //true if we have seen point with same coordinates before
        for(int j=0; j<data.size(); j++){
            if(j==i){continue;} //don't add an edge from data[i] to itself
            Q.push({j,distances[i][j]});
            if(Q.size()>k){
                std::pair<int,double> rem= Q.top();
                Q.pop();
                //keep popped element if it has same distance as top element
                if(rem.second==Q.top().second){
                    Q.push(rem);
                }
            }
            if(distances[i][j]==0 && j<i){
                duplicate=true;
                break;
                //in this case Q has data[j] as first element
            }   
        }
        //now Q contains candidates for data points to be among k closest to data[i], more if those are not unique
        //may have added to many points in attempt to keep points at same distance
        //if duplicate==true we only add edge of weight 0 to first copy for connectedness, no other edges

        if(duplicate){
            std::pair<int,double> nb = Q.top();
            Nodes[nb.first]->mult_++;
            Nodes[i]->mult_=0;
            Nodes[i]->neighbors_[Nodes[nb.first]]=0;
        }
        else{
            //if no duplicate
            //store those not certainly among k closes in uncertainCandidates for a second
            std::vector<std::pair<int,double>> uncertainCandidates;
            int count = Q.size();
            while(count>k){
                std::pair<int,double> nb = Q.top();
                uncertainCandidates.push_back(nb);
                Q.pop();
                count--;
            }
            
            //now Q contains k closest elements, if not unique, 
            //we need to scan uncertainCandidates for elements within  distance Q.top().second;
            if(uncertainCandidates.size()>0){
                double threshold=Q.top().second;
                for(int j = uncertainCandidates.size()-1; j>-1; j--){
                    if(uncertainCandidates[j].second>threshold){break;}
                    Nodes[i]->neighbors_[Nodes[uncertainCandidates[j].first]]=uncertainCandidates[j].second;
                }
            }

               
           
            for( ;!Q.empty();Q.pop()){
                std::pair<int,double> nb = Q.top();
                Nodes[i]->neighbors_[Nodes[nb.first]]=nb.second; //only add this direction, not symmetrically
            }
          
        }
    }
    
    
    Graph* kNN = new Graph(Nodes); 
    return kNN;
}


//adjust for change in mult_ everywhere below


//returns directed kNN graph of G, i.e for each vertex we select the k closest of its direct neighbors
//only selected new neighbors among old neighbors
Graph* getPrekNN(Graph* G, int k){
    //can't take const Graph* G  because of current form of copy constructor for graph
   
    Graph* kNN = new Graph(*G); //kNN is a copy of G
    
    for (int i = 0; i < kNN->size(); i++){   
        //find k nodes closest to currNode below
        Node* currNode = kNN->getVertex(i);
    
        //if currNode has no more than k neighbors anyways, we can just copy all there is nothing to do
        if(currNode->neighbors_.size()<=k){continue;}

        auto cmp = [](std::pair<Node*,double> a, std::pair<Node*,double> b){
            return (a.second<b.second); //this makes Q a max heap
        };
        std::priority_queue<std::pair<Node*,double>, std::vector<std::pair<Node*,double>>, decltype(cmp)> Q(cmp);
        //priority queue storing k closest neighbors

        for(const auto& [key,dist]: currNode->neighbors_){
            Q.push({key,dist});
            if(Q.size()>k){
                std::pair<Node*,double> rem=Q.top();
                Q.pop();
                if(rem.second==Q.top().second){
                    Q.push(rem);
                }
            }
        }
        //now Q contains k closest nodes to currNode
        currNode->neighbors_.clear(); //delete all neighbors, k closest back below
        for (;!Q.empty(); Q.pop()){
            std::pair<Node*,double> temp=Q.top();
            currNode->neighbors_[temp.first]=temp.second;
        }
    }
    return kNN;
}


//helper function to get minimal kNNs, don't call directly!
Graph* getMinkNNHelper(Graph* prekNN){
    //this knows that input graph is a prekNN
    //if nodes in prekNN come with coordinates those are maintained
    //modifies prekNN
    for (int i = 0; i < prekNN->size(); i++)
    {
        Node* currNode=prekNN->getVertex(i);
        auto it = currNode->neighbors_.begin(); //type is {Node*,double}
        while(it != currNode->neighbors_.end()){
            //if currNode not also a neighbor of it->first we delete it->first from neighbors of currNode
            if(it->first->neighbors_.find(currNode)==it->first->neighbors_.end()){
                it=currNode->neighbors_.erase(it);
            }
            else{it++;}
        }    
    }
    return prekNN;
}


//returns minimal kNN, 
//there is an edge between a and b iff a is among k closest neighbors to b and vice versa
Graph* getMinkNN(const std::vector<std::vector<double>>& data, int k){
    return getMinkNNHelper(getPrekNN(data,k));
}

//returns minimal kNN, 
//there is an edge between a and b iff a is among k closest neighbors to b and vice versa
Graph* getMinkNN(Graph* G, int k){
    return getMinkNNHelper(getPrekNN(G,k));
}


//helper function to get maximal kNNs, don't call directly!
Graph* getMaxkNNHelper(Graph* prekNN){
    //this knows that input graph is a prekNN
    //if nodes in prekNN come with coordinates those are maintained
    //modifies prekNN
    for(int i; i<prekNN->size(); i++){
        Node* currNode = prekNN->getVertex(i);
        std::vector<std::pair<Node*,double>> newNeighbors;
        for(auto& [key,dist] : currNode->neighbors_){
            newNeighbors.push_back({key,dist});
        }  
        for(auto [nb, dist]: newNeighbors){
            nb->neighbors_[currNode]=dist;

        }
        /*OLD
        auto it = currNode->neighbors_.begin();
        while(it!=currNode->neighbors_.end()){    
            newNeighbors.push_back(*it);//{it->first,it->second});
            it++;
        }
        {
            it->first->neighbors_[currNode]=it->second; //changes size of it
            it++;
        }
        */
        /*
        for(auto& [key,dist] : currNode->neighbors_){
            key->neighbors_[currNode]=dist;
        }  
        */  
    }
    return prekNN;
}

//returns maximal kNN
//there is an edge between a and b if a is among k closest neighbors to a or b among k closest nb to a
Graph* getMaxkNN(Graph* G, int k){
    return getMaxkNNHelper(getPrekNN(G,k));
}

//returns maximal kNN
//there is an edge between a and b if a is among k closest neighbors to a or b among k closest nb to a
Graph* getMaxkNN(const std::vector<std::vector<double>>& data, int k){
    Graph* prekNN = getPrekNN(data,k); //on heap.
    Graph* maxkNN=getMaxkNNHelper(prekNN); //modifies graph above
    return maxkNN;
}