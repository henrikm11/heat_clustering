//graph.cpp
//contains implementation of some graph algorithms tailored towards our application
//graphs undirected? but weighted?

/*
TODO:
1D clustering?
figuring out time scales for heat flow - compare to diameter etc.
decomposing into connected components
*/

#include "clustering.h"



//constructors for Node
Node::Node(void): mult_(1){return;}
Node::Node(int mult): mult_(mult){return;}
Node::Node(int d, int mult): d_(d), mult_(mult){return;};

Node::Node(const std::vector<double>& coordinates, int mult): coordinates_(coordinates), d_(coordinates.size()), mult_(mult){return;}
Node::Node(const std::vector<double>& coordinates, const std::unordered_map<Node*,double>& neighbors, int mult):coordinates_(coordinates), d_(coordinates.size()), neighbors_(neighbors), mult_(mult){return;}
Node::Node(const std::unordered_map<Node*,double>& neighbors, int mult):neighbors_(neighbors),mult_(mult){return;} // use any neighbor if existent to set d_>


//constructors for Graph
Graph::Graph(void):labeled_(false),distanceMatrixInitialized_(false), distancesComputed_(false), connectedChecked_(false){return;}
Graph::Graph(std::vector<Node*> vertices): vertices_(vertices),labeled_(false),distanceMatrixInitialized_(false), distancesComputed_(false), connectedChecked_(false){return;}

//would like to change to const Graph&, but we may need to create labels if necessary
Graph::Graph(Graph& G)
    :labeled_(G.labeled_), distanceMatrixInitialized_(G.distanceMatrixInitialized_), distancesComputed_(G.distancesComputed_),distances_(G.distances_),
    connected_(G.connected_),connectedChecked_(G.connectedChecked_),connectedCompCount_(G.connectedCompCount_)
    {
        for(int i=0; i<G.size(); i++){
            Node* temp = new Node;
            this->vertices_.push_back(temp);
            this->label_[temp]=i;
        }
        this->labeled_=true;
        //have created Nodes with labels, now assign neighbors
        G.labelNodes(); //make sure G also has labels
        for(int i=0; i<G.size(); i++){
            Node* currNode=this->getVertex(i);
            Node* currNodeOld=G.getVertex(i);
            for(const auto& edge : currNodeOld->neighbors_){
                int nbLabel=G.getLabel(edge.first);
                currNode->neighbors_[this->vertices_[nbLabel]]=edge.second;
            }
        }
        //have copied all edges
        if(G.connectedChecked_ && !(G.connected_)){
            for(const auto& [node,value] : G.componentLabel_){
                int idx=G.getLabel(node);
                this->componentLabel_[this->vertices_[idx]]=value;
            }
        }
        return;
    }

Graph& Graph::operator=(Graph& G){
    this->labeled_=G.labeled_;
    this->distanceMatrixInitialized_=G.distanceMatrixInitialized_;
    this->distancesComputed_=G.distancesComputed_;
    this->distances_=G.distances_;
    this->connected_=G.connected_;
    this->connectedChecked_=G.connectedChecked_;
    this->connectedCompCount_=G.connectedCompCount_;
    this->vertices_={};
    for(int i=0; i<G.size(); i++){
        Node* temp = new Node;
        this->vertices_.push_back(temp);
        this->label_[temp]=i;
        }
    for(int i=0; i<G.size(); i++){
            Node* currNode=this->getVertex(i);
            Node* currNodeOld=G.getVertex(i);
            for(const auto& edge : currNodeOld->neighbors_){
                int nbLabel=G.getLabel(edge.first);
                currNode->neighbors_[this->vertices_[nbLabel]]=edge.second;
            }
        }
    if(G.connectedChecked_ && !(G.connected_)){
        for(const auto& [node,value] : G.componentLabel_){
            int idx=G.getLabel(node);
            this->componentLabel_[this->vertices_[idx]]=value;
        }
    }
    return *this;
}


Graph::~Graph(void){} //needs to be modified later on


//basic member variable access functions for Graph
int Graph::size()const{return this->vertices_.size();}

Node* Graph::getVertex(int i) const{
    if (i<0 || i>=this->size()){
        throw std::out_of_range("Graph::getVertex()");
    }
    return this->vertices_[i];
}

std::vector<Node*>* Graph::getVertices(){
    return &(this->vertices_);
}

double Graph::getDistance(Node* source, Node* target){
    //assume graph is connected for now
    //this is running Dijkstra if necessary
    double inf = std::numeric_limits<double>::infinity();
    // std::unordered_map<Node*, int> labels=this->getLabel();
    //make sure we have associated lables
    if(!distanceMatrixInitialized_){    
        std::vector<double> temp(this->vertices_.size(),inf);
        this->distances_= std::vector<std::vector<double>>(this->vertices_.size(), temp);
        distanceMatrixInitialized_=true;
    }
    else{
        //can check if we have done this computation before
        double dist=this->distances_[this->getLabel(source)][this->getLabel(target)];
        if(dist<inf){
            return dist;
        }
    }

    //now we know that we have labels and the matrix distances_  

    auto cmp = [] (std::pair<int,double> a, std::pair<int,double> b) {return (a.second>b.second);};   
    std::priority_queue<std::pair<int,double>, std::vector<std::pair<int,double>>,decltype(cmp)> Q(cmp);  
    //queue of vertices still to be checked
    //stored as (int label, double distance), compared by distance
    //use lazy removal for updates

    //push neighbors of source onto Q
    //initialized with distance by edge connecting them
    for(const auto& [key, value] : source->neighbors_){
        std::pair<int,double>edge={this->getLabel(key),value};
        Q.push(edge);
    }

    std::unordered_map<int, bool> checked; 
    //true at i iff distance has been computed to vertices_[i]
    int sourceLabel;
    for (int i = 0; i < this->vertices_.size(); i++){
        if(vertices_[i]==source){
            checked[i]=true;
            sourceLabel=i;
        }
        checked[i]=false;
    }
    
    int targetLabel=this->getLabel(target);  
    while(!checked[targetLabel]){
        //find vertex that we can check off in this round
        std::pair <int,double> checkNext = Q.top();
        while(checked[checkNext.first]==true){
            //because of lazy removal, there may be vertices in the queue that have already been checked
            Q.pop();
            checkNext=Q.top();
        }
        //now checkNext is the vertex we can check in this step
        this->distances_[checkNext.first][sourceLabel]=checkNext.second;
        this->distances_[sourceLabel][checkNext.first]=checkNext.second;
        checked[checkNext.first]=true;
        

        //push neighbors of checkNext that have not been checked yet onto Q
        for(const auto& [key, value] : this->vertices_[checkNext.first]->neighbors_){
            if( !checked [this->getLabel(key) ] ){
                std::pair<int,double>edge={this->getLabel(key),checkNext.second+value};
                //push with distance being distance(source, checkNext)+this edge
                Q.push(edge);
            } 
        }
    }

    return this->distances_[sourceLabel][targetLabel];
}

int Graph::getLabel(Node* a){
    if(!this->labeled_){
        this->labelNodes();
    }
    return this->label_[a];
}

/*
std::unordered_map<Node*, int> Graph::getLabel(void){
    if(!this->labeled_){
        this->labelNodes(); 
    }
    return this->label_;
}
*/

bool Graph::isConnected(){
    this->checkConnected_();
    return this->connected_;
}

int Graph::connectedCompCount(){
    if(this->isConnected()){return 1;}
    this->labelComponents_();
    return this->connectedCompCount_;
}

int Graph::connectedCompLabel(Node* a){
    this->labelComponents_(); //make sure we are not doing this over and over again
    int lab=this->componentLabel_[a];
    return lab;
}



//basic functions changing private member variables of Graph

void Graph::insert(Node* node){
    //should we make sure that this a new node?
    this->vertices_.push_back(node);
    if(this->labeled_){
        this->label_[node]=this->label_.size();
    }
    return;
}

void Graph::insert(Node* a,Node* b, double weight){
    a->neighbors_[b]=weight;
    b->neighbors_[a]=weight;
    return;
}

void Graph::remove(Node* a, Node* b){
    if(b->neighbors_.count(a)==0){
        return;
    }
    b->neighbors_.erase(a);
    a->neighbors_.erase(b);
    return;
} 

void Graph::labelNodes(){
    int s= this->size();
    if(this->label_.size()==s){return;} 
    for (int i = 0; i < s; i++)
    {
        this->label_[this->vertices_[i]]=i;
    }
    return;
}

void Graph::computeDistances(){
    //assume graph is connected for now
    //this is running Dijkstra if necessary
    //there must be a way of doing this more efficiently using computations one has already done
    //improvement should be by factor 2
    if(this->distancesComputed_){return;}
    if(!distanceMatrixInitialized_){   
        double inf = std::numeric_limits<double>::infinity();
        // std::unordered_map<Node*, int> labels=this->getLabel();
        std::vector<double> temp(this->vertices_.size(),inf);
        this->distances_= std::vector<std::vector<double>>(this->vertices_.size(), temp);
        distanceMatrixInitialized_=true;
    }
    //now we know that we have labels and the matrix distances_  

    int s = this->size();
    for(int i=0; i<s; i++){
        //compute all distances from source vertices_[i]
        auto cmp = [] (std::pair<int,double> a, std::pair<int,double> b) {return (a.second>b.second);};   
        std::priority_queue<std::pair<int,double>, std::vector<std::pair<int,double>>,decltype(cmp)> Q(cmp);  
        //queue of vertices still to be checked
        //stored as (int label, double distance), compared by distance
        //use lazy removal for updates

        Node* source = this->vertices_[i];
        //push neighbors of source onto Q
        //initialized with distance by edge connecting them
        for(const auto& [key, value] : source->neighbors_){
            std::pair<int,double>edge={this->getLabel(key),value};
            Q.push(edge);
        }
        std::unordered_map<int, bool> checked; 
        //true at i iff distance has been computed to vertices_[i]
        int sourceLabel;
        for (int i = 0; i < this->vertices_.size(); i++){
            if(vertices_[i]==source){
                checked[i]=true;
                sourceLabel=i;
            }
            checked[i]=false;
        }

        int count=1; //counts number of vertices for which we have computed distance already
        while (count<s){
            //find vertex that we can check off in this round
            std::pair <int,double> checkNext = Q.top();
            while(checked[checkNext.first]==true){
                //because of lazy removal, there may be vertices in the queue that have already been checked
                Q.pop();
                checkNext=Q.top();
            }
            //now checkNext is the vertex we can check in this step
            this->distances_[checkNext.first][sourceLabel]=checkNext.second;
            this->distances_[sourceLabel][checkNext.first]=checkNext.second;
            checked[checkNext.first]=true;
            count++;

            //push neighbors of checkNext that have not been checked yet onto Q
            for(const auto& [key, value] : this->vertices_[checkNext.first]->neighbors_){
                if( !checked [this->getLabel(key) ] ){
                    std::pair<int,double>edge={this->getLabel(key),checkNext.second+value};
                    //push with distance being distance(source, checkNext)+this edge
                    Q.push(edge);
                } 
            }
        }   
    }

    this->distancesComputed_=true;
    return;
}


void Graph::helperDFS_(
    Node* source,
    int currComponentLabel,
    std::unordered_map<Node*,bool>& checkedDFS,
    int& count
    )
    {
    checkedDFS[source]=true;
    count++;
    this->componentLabel_[source]=currComponentLabel;
    for(auto& [key,value]: source->neighbors_){
        if(checkedDFS.count(key)==0){
            helperDFS_(key, currComponentLabel, checkedDFS, count);
        }
    }
    return;
}

void Graph::checkConnected_(){
    if(this->connectedChecked_){return;}
    if(this->size()==0){
        this->connectedChecked_=true;
        this->connected_=true;
        return;
    }

    Node* source = this->getVertex(0);
    std::unordered_map<Node*,bool> checkedDFS;
    int count=0;

    this->helperDFS_(source, 0, checkedDFS, count);

    if(count==this->size()){
        this->connected_=true;
    }
    else{
        this->connected_=false;
    }

    this->connectedChecked_=true;
    return;
}


void Graph::labelComponents_(){
    this->checkConnected_();
    if(this->connected_){
        this->connectedCompCount_=1;
        return;
    }

    int currComponentLabel = 0;
    std::unordered_map<Node*,bool> checkedDFS;
    int count = 0;

    for (int i = 0; i < this->size(); i++)
    {
        if(checkedDFS.count(this->getVertex(i))==0){
            //haven't been here before
            this->helperDFS_(this->getVertex(i), currComponentLabel, checkedDFS,count);
            currComponentLabel++;
        }   
    }
    this->connectedCompCount_=currComponentLabel; //+1?
    return;  
}




