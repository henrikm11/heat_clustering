//node.cpp
//contains implementation of node struct

#include "clustering.h"
#include "graph.h"

Node::Node()
    :neighbors_(0),
    mult_(1),
    d_(-1),
    coordinates_(0)
    {}

Node::Node(const int mult)
    :neighbors_(0),
    mult_(mult),
    d_(-1),
    coordinates_(0)
    {}

Node::Node(const int d, int mult)
    :neighbors_(0),
    mult_(mult),
    d_(d),
    coordinates_(0)
    {}

Node::Node(const std::vector<double>& coordinates, int mult)
    :neighbors_(0),
    mult_(mult),
    d_(coordinates.size()),
    coordinates_(coordinates)
    {}

Node::Node(const std::vector<std::pair<Node*,double>>& neighbors, int mult)
    :neighbors_(neighbors),
    mult_(mult),
    d_(-1),
    coordinates_(0)
    {} 

Node::Node(const std::vector<double>& coordinates,const std::vector<std::pair<Node*,double>>& neighbors, int mult)
    :neighbors_(neighbors),
    mult_(mult),
    d_(coordinates.size()),
    coordinates_(coordinates)
    {}


/*
Node::Node(const Node& originalNode)
    :neighbors_(originalNode.neighbors_),
    mult_(originalNode.mult_),
    d_(originalNode.d_),
    coordinates_(originalNode.coordinates_)
    {}
*/



/*
Node::~Node(){
    
    std::cout << this->coordinates_[0] <<",";
    std::cout << this->neighbors_.size()<<std::endl;
    for(auto& it1 : neighbors_){
        //std::cout <<"here" << std::endl;;
        Node* nb = it1.first;
        double dist = it1.second;
        if(nb==nullptr) continue;
        //if(nb!=nullptr){
        //make sure this hasn't been invalidated before 
        std::pair<Node*,double> edge = {this, it1.second};
        int s=nb->neighbors_.size();
        auto it2=nb->neighbors_.begin();
        if(it2!=nb->neighbors_.end()){
            nb->neighbors_.erase(it2);
        }
        
        for(int pos=0; pos<s; pos++){
            if(*it2==edge){
                //it2=nb->neighbors_.erase(it2);
                break;
            }
            it2++;
        }
        for(auto it2 = nb->neighbors_.begin(); it2!=nb->neighbors_.end();){
            if(*it2==edge){
                std::cout <<"here" << std::endl;;
                it2=nb->neighbors_.erase(it2);
                break;
            }
            else{
                it2++;
            } 
        }   
            
        //}
    }
 
}
*/




std::vector<double> Node::getCoordinates(void) const{
    return coordinates_;
}

void Node::addNeighbor(Node& nb, double dist){
    Node* nbPtr =&nb;
    std::pair<Node*,double> edge{nbPtr,dist};
    addNeighbor(edge);
    return;
}
/*
    addNeighbor(<)
    //check if neighbor already present
    for(size_t i=0; i<neighbors_.size();i++){
        if(neighbors_[i].first==&nb){
            neighbors_[i].second=dist;
            return;
        }
    }
    Node* nbPtr=&nb;
    std::pair<Node*,double> edge{nbPtr,dist};
    neighbors_.push_back(edge);
    return;
    
}
*/

 void Node::addNeighbor(std::pair<Node*,double> edge){
    if(edge.first==nullptr) return;
    Node* nbPtr=edge.first;
    double dist=edge.second;
  
    for(size_t i=0; i<neighbors_.size();i++){
        if(neighbors_[i].first==nbPtr){
            neighbors_[i].second=dist;
            return;
        }
    }
    neighbors_.push_back(edge);
    return;
 }


//return pointer arithmetic is nonsense
std::vector<std::pair<Node*,double>>::iterator Node::removeNeighborIt(Node& removeNode){
    Node* remNodePtr=&removeNode;
    auto it=neighbors_.begin();
    for(auto it=neighbors_.begin(); it!=neighbors_.end(); ){
        Node* currNb=it->first;
        if(remNodePtr==currNb){
            it=neighbors_.erase(it);
            break;
        }
        it++;
    }
    return it;
}


void Node::removeNeighbor(Node& node){
    
    Node* nodePtr=&node;
    auto it = std::remove_if(
        neighbors_.begin(), 
        neighbors_.end(),
        [nodePtr](const std::pair<Node*,double>& edge){return edge.first==nodePtr;}
    );
    
    neighbors_.erase(it, neighbors_.end());
    return;   
 }

 bool Node::checkNeighbor(const Node& node) const{
    for(const auto& edge : neighbors_){
        if(edge.first==&node) return true;
    }
    return false;
 }