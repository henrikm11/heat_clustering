//heatFlow.cpp
//contains implementation of simulation of heat flow on graphs

#include "clustering.h"

//compute graph laplacian of given functions
std::unordered_map<Node*, double> Graph::laplacian(const std::unordered_map<Node*, double>& func, bool& admissible){
    //built in check if size of vertices of graph matches size of func
    int s = this->size();
    int count=0;
    if(!admissible)
    {
        //check if function is defined on given graph
        std::unordered_map<Node*,int> vertexCount;
        for (int i = 0; i < s; i++)
        {
            if(vertexCount.count(this->vertices_[i])>0 || count > s)
            {
                std::cout << "Function is not defined on graph" << std::endl;
                exit(1);
            }
            vertexCount[this->vertices_[i]]=true;
            count+=1;
        }
    }
    //if we make it to here function is defined on actual graph, no need to handle this again in future operations

    admissible=true;
    double comp;
    std::unordered_map<Node*, double> lap;
    for (int i = 0; i < this->size(); i++)
    {    
        Node* currNode=this->getVertex(i);
        lap[this->getVertex(i)]=0;
        for(const auto& [key,value] :currNode->neighbors_)
        {
            //[key,value] is neighbouring Node*, weight
            double con;
            if (value!=0){con=1/value;}
            else{con=0;} //these terms are taken into account through mult_
            lap[this->vertices_[i]]+=
            ((key->mult_)*func.at(key)-(currNode->mult_)*func.at(currNode))*con;
        }
        //this is taking into account weights through term con  
    }
    return lap;
}


//computes on time discrete step in Euler iteration to solve heat equation
//evolves heat distribution given by func in direction laplacian(func) by timeScale
void Graph::heatIterationStep(double timeScale, std::unordered_map<Node*,double>& func, bool& admissible){
    int s = this->size();
    std::unordered_map<Node*, double> lap=this->laplacian(func, admissible);
    for (int i = 0; i < s; i++)
    {
        func[this->vertices_[i]]+=(timeScale*lap[this->vertices_[i]]);
    }
    return;
}

//computes heat evolution over time using approximations of size timeScale
void Graph::heatDiffusion(double time, double timeScale, std::unordered_map<Node*,double>& func, bool& admissible){
    int k=time/timeScale+1;
    for (int i = 0; i < k; i++)
    {
        this->heatIterationStep(timeScale, func, admissible);
    }
    return;
}


