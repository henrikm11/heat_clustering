//heatflow.cpp

#include "clustering.h"
#include "graph.h"

//TODO
/*
dynamic time steps for iterations to have max principle
this is sufficient to ensure convergence to equilibrium
*/


std::vector<double> Graph::laplacian(const std::vector<double>& func){
    if(func.size()!=size()){
        //func not a well defined function on graph
        throw(std::domain_error("Graph::laplacian"));
    }
    
    std::vector<double> lap(size(),0);  
    for(size_t i =0; i<size(); i++){
        //compute laplacian in vertices_[i]
        int mult = vertices_[i]->mult_;
        for(size_t j =0; j<adjacencyLists_[i].size();j++){
            int ind = adjacencyLists_[i][j].first;
            double dist = adjacencyLists_[i][j].second;
            double weight=0;
            int multNb=vertices_[i]->mult_;
            if(dist!=0){weight=1/dist;}
            lap[i]+=(multNb*func.at(ind)-mult*func.at(i))*weight; //check formula!
        }
    }  

    return lap;
}

void Graph::heatIterationStep(std::vector<double>& initialDist, double timeScale, double& low, double& high){
    std::vector<double> lap=laplacian(initialDist);

    if(low==-1){low=*(std::min_element(initialDist.begin(),initialDist.end()));}
    if(high==-1){high=*(std::max_element(initialDist.begin(),initialDist.end()));}

    bool maxPrinciple=true;
    double dynamicTimeScale=timeScale;
    double newHigh=low; //stores max after update
    double newLow=high; //stores min after update

    

    for(size_t i =0; i<size(); i++){
        initialDist[i]+=timeScale*lap[i];
        newHigh=std::max(newHigh, initialDist[i]);
        newLow=std::min(newLow, initialDist[i]);
        if(initialDist[i]>high){
            maxPrinciple=false;
            double beforeUpdate=initialDist[i]-timeScale*lap[i];
            dynamicTimeScale=std::min(dynamicTimeScale, (high-beforeUpdate)/(2*lap[i]));
        }
        else if(initialDist[i]<low){
            maxPrinciple=false;
            double beforeUpdate= initialDist[i]-timeScale*lap[i];
            dynamicTimeScale=std::min(dynamicTimeScale, (low-beforeUpdate)/(2*lap[i]));
        }
    }
    if(!maxPrinciple){
        //undo iteration
        for(size_t i=0; i<size(); i++){
            initialDist[i]-=timeScale*lap[i];
        }
        //run iteration with adjusted smaller timeScale
        heatIterationStep(initialDist, dynamicTimeScale, low, high);
        //this also updates low and high accordingly
        
    }
    else{
        //everything went fine up to here
        low=newLow;
        high=newHigh;
    }
    return;
}