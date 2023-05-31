//OneDimClustering.cpp

#include "clustering.h"

//1Dclustering.cpp
//implementation of clustering in 1 dimension
//specifically implements the function oneDimensionalClustering, everything else is helper functions
//we apply this to Graph after dimension reduction via heat flow


double OneDimClusterHelper::getBaseDensity(const std::vector<std::pair<int,double>>& dataInd){
    if(dataInd.size()==0){return -1;}
    if(dataInd.size()<100){
        if(dataInd.front().second==dataInd.back().second){return -1;}
        double baseDensity=1/(dataInd.back().second-dataInd.front().second);
    }
    int start=0;
    int end=dataInd.size()-1;
    while(start<end && (start+1)<0.01*dataInd.size()){
        start++;
    }
    while(start<end && dataInd.size()-end<0.01*dataInd.size()){
        end--;
    }
    //now 0.99 of data are to the right of start
    //0.99 of data are to the left of end
    if(dataInd[start].second==dataInd[end].second){return -1;}
    double baseDensity=1/(dataInd[end].second-dataInd[start].second);
    return baseDensity;
}

double OneDimClusterHelper::getExpectedBaseWindowCount(const std::vector<std::pair<int,double>>& dataInd, const double bandWidth){
    double baseDensity=getBaseDensity(dataInd);
    if(baseDensity==-1){return dataInd.size();} //data essentially concentrated at one point
    double expectedCount=bandWidth*dataInd.size()*baseDensity;
    return expectedCount;
}

std::vector<std::pair<int,int>> OneDimClusterHelper::getClusterSplits(const std::vector<std::pair<int,double>>& dataInd, const double significance, const double bandWidth, const double minClusterSize){
    assert(significance>0);
    assert(bandWidth>0);
    assert(minClusterSize>0 && minClusterSize<=1);
    
    
    std::vector<std::pair<int,int>> clusterSplits;
    clusterSplits.clear(); //

    //iteration parameters for sliding window
    int centerOfWindow=dataInd.size()-1;
    int endOfWindow=centerOfWindow;
    int startOfWindow=centerOfWindow;
    double dataCenter=dataInd[dataInd.size()-1].second;
    int steps=0;
    int lastClusterSplit=dataInd.size()-1;

    //scan for places in which density is below significance*maxDensity
    //keep track if point with high density in between
    int currLabel=0;
    bool cluster=false;
    int count=0;

    double expectedBaseWindowCount = getExpectedBaseWindowCount(dataInd, bandWidth);
    if(expectedBaseWindowCount==dataInd.size()){
        clusterSplits.push_back({0,dataInd.size()-1});
        return clusterSplits;
    }
    

    double inClusterThreshold = significance*expectedBaseWindowCount;
    double splitClusterThreshold =expectedBaseWindowCount/significance;

    while(centerOfWindow>-1)
    {
        if(steps%2==0)
        {
            dataCenter=dataInd[centerOfWindow].second;
        }
        else
        {
            dataCenter=(dataInd[centerOfWindow+1].second+dataInd[centerOfWindow].second)/2;
        }

        //compute density in window of size bandWidth centered at centerOfWindow
        while(dataInd[endOfWindow].second-dataCenter>bandWidth)
        {
            endOfWindow--;
        }
        //now window ends at endOfWindow (inclusive)
        while(startOfWindow>0 && dataCenter-dataInd[startOfWindow-1].second<=bandWidth)
        {
            startOfWindow--;
        }
        //now window starts at startOfWindow (inclusive)
        int currDensity=endOfWindow-startOfWindow+1;
        if(!cluster && currDensity>=inClusterThreshold){cluster=true;}
        if(cluster && currDensity<=splitClusterThreshold && lastClusterSplit-centerOfWindow+1>=minClusterSize*dataInd.size())
        {
            //detected cluster     
            std::pair<size_t,size_t> currSplit{centerOfWindow, lastClusterSplit};
            clusterSplits.push_back(currSplit);
            lastClusterSplit=centerOfWindow-1;
            currLabel++;
            cluster=false;
        }
        //update of iteration parameters
        if(steps%2==0)
        {
            centerOfWindow--;
        }
        steps++;
    }
    if(cluster){
        clusterSplits.push_back({0,lastClusterSplit});
    }
    return clusterSplits;
}



std::vector<int> oneDimensionalClustering(const std::vector<double>& data, const double significance, const double bandWidth, const double minClusterSize){
    if(data.size()==0){return std::vector<int>(0);}
    if(significance<0){
        throw(std::domain_error("oneDimensionalClustering: significance"));
    }
    
    std::vector<std::pair<int,double>> dataInd; //{i,data[i]}
    dataInd.reserve(data.size());
    for(size_t i=0; i<data.size(); i++){
        dataInd.push_back({i,data[i]});
    }
    //sort data
    auto cmp = [] (const std::pair<int,double>& a, const std::pair<int,double>& b){return a.second<b.second;};
    std::sort(dataInd.begin(),
        dataInd.end(),
        cmp);
    
    //Helper object
    OneDimClusterHelper Helper;

    //get indices at which clusters begin/end
    std::vector<std::pair<int,int>> clusterSplits = Helper.getClusterSplits(
        dataInd,
        significance,
        bandWidth, 
        minClusterSize
        );
    
    //label data according to cluster splits
    std::vector<int> clusterLabels(dataInd.size(),-1);
    int currLabel=0;
    for(const auto& split : clusterSplits){
        size_t start=split.first;
        size_t end=split.second;
        for(size_t i=start; i<=end; i++){
            clusterLabels[dataInd[i].first]=currLabel;
        }
        currLabel++;
    }
    return clusterLabels;
}