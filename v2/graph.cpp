//graph.cpp
//contains implementation of
//class graph

#include "clustering.h"
#include "graph.h"

//constructors

Graph::Graph()
    :vertices_(0),
    adjacencyLists_(0)
    {}

Graph::Graph(const std::vector<Node*>& vertices)
    :vertices_(vertices.size()),
    adjacencyLists_(vertices_.size())
    {
    
        //create new Nodes on heap
        std::unordered_map<Node*,int> originalLabels;
        for(size_t i=0; i<vertices.size(); i++){
            Node* originalNode = vertices[i];
            originalLabels[originalNode]=i;
            vertices_[i] = new Node(*originalNode); //this has neighbors in vertices not vertices_!
            vertices_[i]->neighbors_.clear();
        }
        //fix neighbors
        //also checks if any neighbors not part of graph, if not we don't insert
        for(size_t i=0; i<vertices_.size(); i++){
            Node* originalNode = vertices[i];
            Node* newNode = vertices_[i];
            for(const auto& edge : originalNode->neighbors_){
                if(originalLabels.count(edge.first)==0){continue;} //neighbor not contained in graph 
                Node* newNb = vertices_[originalLabels[edge.first]];
                newNode->addNeighbor({newNb, edge.second});
                newNb->addNeighbor({newNode, edge.second});
            }
        } 

        std::unordered_map<Node*,int> newLabels;
        for(size_t i=0; i<vertices_.size(); i++){
            newLabels[vertices_[i]]=i;
        }
        for(size_t i=0; i<vertices_.size(); i++){
            for(const auto& nb : vertices_[i]->neighbors_){
                adjacencyLists_[i].push_back({newLabels[nb.first],nb.second});
            }
        }

    }
//end of constructor


//copy and assignment

Graph::Graph(const Graph& originalGraph)
    :Graph::Graph(originalGraph.getVertices())
    {}

Graph& Graph::operator=(const Graph& originalGraph){
    //delete old nodes
    for(size_t i=0; i<size(); i++){
        delete vertices_[i];
        vertices_[i]=nullptr;
    }
    vertices_.clear();
    vertices_.reserve(originalGraph.size());

    std::unordered_map<Node*,int> originalLabels;
    for(size_t i=0; i<originalGraph.size(); i++){
        Node* originalNode = &(originalGraph.getVertex(i));
        originalLabels[originalNode]=i;
        vertices_.push_back(new Node(*originalNode)); //this has neighbors in original not copy
    }

    //fix neighbors
    //by hand more efficient than calling addNeighbor since originalGraph had no redundancies
    for(size_t i=0; i<originalGraph.size(); i++){
        Node* originalNode = &(originalGraph.getVertex(i));
        Node* newNode = vertices_[i];
        for(const auto& [nb,dist] : originalNode->neighbors_){
            Node* newNb = vertices_[originalLabels[nb]];
            newNode->neighbors_.push_back({newNb,dist});
        }
    } 

    adjacencyLists_=originalGraph.adjacencyLists_;

    return *this;
}

//destructor
Graph::~Graph(){
    for(size_t i=0; i<size(); i++){
        delete vertices_[i];
        vertices_[i]=nullptr;
    }
}



//basic member access functions
int Graph::size()const{
    return vertices_.size();
}

Node& Graph::getVertex(int i) const {
    if (i<0 || i>=size())
    {
        throw std::out_of_range("Graph::getVertex()");
    }
    return *(vertices_[i]);
}

const std::vector<Node*>& Graph::getVertices() const{
    return vertices_;
 }


//functions changing members

int Graph::findNode(const Node& node) const {
    int idx=-1; //default return if not presen
    for(size_t i =0; i<size(); i++){
        if(vertices_.at(i)==&node){
            idx=i;
            break;
        }
    }
    return idx;
}

void Graph::insert(const Node& newNode){

    //check if newNode already present
    if(findNode(newNode)!=-1) return;

    //insert copy of newNode
    Node* insertedNode = new Node(newNode);
    std::vector<std::pair<Node*,double>> potentialNeighbors=insertedNode->neighbors_;
    insertedNode->neighbors_.clear();
    vertices_.push_back(insertedNode);

    //add edges
    std::unordered_map<Node*,int> reqLabels;
    for(const auto& edge : potentialNeighbors){
        //check if end of edge is contained in Graph
        int edgeIdx = findNode(*edge.first);
        if(edgeIdx==-1){std::cout<<"ok";continue;} //end of edge not in graph

        assert(edgeIdx>=0);
        assert(edgeIdx<size());
        
        reqLabels[vertices_[edgeIdx]]=edgeIdx;
        edge.first->addNeighbor({insertedNode,edge.second});
        insertedNode->addNeighbor({edge.first,edge.second});
    }
    //update adjacencyLists_
    int newIdx = vertices_.size()-1; //index of newly added node
    std::vector<std::pair<int,double>> newAdj; //adjacency for newly added node
    for(const auto& edge : insertedNode->neighbors_){
        int idx = reqLabels[edge.first];
        adjacencyLists_[idx].push_back({newIdx, edge.second});
        newAdj.push_back({idx, edge.second});   
    }
    adjacencyLists_.push_back(newAdj);

    return;
}

bool Graph::checkNeighbors(const Node& node1, const Node& node2, bool NodesOK )const {
    if(!NodesOK){
        //nodes may not even be in graph
        if(findNode(node1)==-1) return false;
        if(findNode(node2)==-1) return false;
    }
    
    if(node1.neighbors_.size()<node2.neighbors_.size()){
        return node1.checkNeighbor(node2);
    }
    return node2.checkNeighbor(node1);
}

int Graph::countNeighbors(const Node& node) const {
    return node.neighbors_.size();
}

int Graph::countNeighbors(int i) const{
    if(i<0 || i>=size()){
        throw(std::out_of_range("Graph::countNeighbors"));
    }
    return vertices_[i]->neighbors_.size();
}

void Graph::insert(Node& node1, Node& node2, double dist){
    int idx1 = findNode(node1);
    int idx2 = findNode(node2);
    //exit if one of the nodes not in graph
    if(idx1==-1 || idx2==-1) return;

    bool nbBefore = checkNeighbors(node1, node2);

    if(nbBefore){
        //only update distance in nodes and adjacency
        for(auto& edge : node1.neighbors_){
            if(edge.first == &node2){
                edge.second = dist;
                break;
            }
        }
        for(auto& edge : node2.neighbors_){
            if(edge.first == &node1){
                edge.second = dist;
                break;
            }
        }
        for(auto& edge : adjacencyLists_[idx1]){
            if(edge.first == idx2){
                edge.second = dist;
                break;
            }
        }
        for(auto& edge : adjacencyLists_[idx2]){
            if(edge.first == idx1){
                edge.second = dist;
                break;
            }
        }
    }
    else{
        node1.addNeighbor(node2, dist);
        node2.addNeighbor(node1,dist);
        adjacencyLists_[idx1].push_back({idx2,dist});
        adjacencyLists_[idx2].push_back({idx1,dist});
    }

    return;
}

void Graph::remove(Node* node){
    if(node==nullptr){return;}
    int idx = findNode(*node);
    if(idx==-1){
        throw(std::domain_error("Graph::remove, attempting to delete node not in graph"));
    }
    for(auto& edge : node->neighbors_){
        //remove node as neighbors from everyhwere
        remove(*node, *edge.first);
    }
    auto it =  std::remove(vertices_.begin(), vertices_.end(), node);
    vertices_.erase(it, vertices_.end());
    delete node;
    return;
}

void Graph::remove(Node& node1, Node& node2){
    //remove neighbors from nodes
    node1.removeNeighbor(node2);
    node2.removeNeighbor(node1);

    //update adjacency list
    int idx1=findNode(node1);
    int idx2=findNode(node2);

    auto it1 = std::remove_if(adjacencyLists_[idx1].begin(),
            adjacencyLists_[idx1].end(),
            [idx2](std::pair<int,double> edge){return edge.first==idx2;}
        );
    adjacencyLists_[idx1].erase(it1, adjacencyLists_[idx1].end());

    auto it2 = std::remove_if(adjacencyLists_[idx2].begin(),
            adjacencyLists_[idx2].end(),
            [idx1](std::pair<int,double> edge){return edge.first==idx1;}
        );
    adjacencyLists_[idx2].erase(it2, adjacencyLists_[idx2].end());
    return;
}


double Graph::getMinDistance(){
    if(size()<=1){return 0;}
    double minDistance=-1;
    for(size_t i=0; i<size();i++){
        for(const auto&  edge: adjacencyLists_[i]){
            if(minDistance==-1){minDistance=edge.second;}
            minDistance=std::min(minDistance, edge.second);
        }
    }
    if(minDistance==-1){minDistance=0;}
    return minDistance;
}

void Graph::normalizeDistances(double newMinDist){
    double minDistance=getMinDistance();
    if(minDistance==0){return;}
    for(size_t i=0; i<size();i++){
        for(auto& edge: adjacencyLists_[i]){
            edge.second/=minDistance;
            edge.second*=newMinDist;
        }
    }
    for(size_t i=0; i<size();i++){
        for(auto& edge : vertices_[i]->neighbors_){
            edge.second/=minDistance;
            edge.second*=newMinDist;
        }
    }
    return;
}

bool Graph::isConnected(){
    if(size()<=1){return true;}
    std::set<Node*> visited;
    std::queue<Node*> visitNext({&getVertex(0)});
    while(!visitNext.empty()){
        Node* currentNode=visitNext.front();
        visited.insert(currentNode);
        visitNext.pop();
        for(const auto& edge : currentNode->neighbors_){
            if(visited.count(edge.first)>0){continue;}
            visitNext.push(edge.first);
        }
    }
    return visited.size()==size();
}

bool Graph::checkSym(){
    for(const auto& v : vertices_){
        for(const auto& w : vertices_){
            bool check1=v->checkNeighbor(*w);
            bool check2=w->checkNeighbor(*v);
            if(check1 && !check2){return false;}
            if(check2 && !check1){return false;}
        }
    }
    return true;

}
