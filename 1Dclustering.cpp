//1Dclustering.cpp
//implementation of clustering in 1 dimension
//we apply this to Graph after dimension reduction via heat flow

#include "clustering.h"


//attempts 1D clustering via kernel density estimator
std::unordered_map<double,int>* oneDimensionalClustering(std::vector<double>& data, double significance, double bandWidth, double minClusterSize){
    //try to locate clusters by locating jumps in density
    //approximates (up to normalization) density of one dimensional distribution using step functions of width bandWidth
    //identifies clusters if points of large density are separated by points of small density
    //if maxDensity>=significance*minDensity we assign to a cluster
    //separation between different clusters if density drops low aain
   
    int currLabel=0;
    std::unordered_map<double, int>* clusterLabel = new std::unordered_map<double, int>;

    std::sort(data.begin(),data.end());
    
    //compute density in windows of size bandwidth centered at every data point
    int centerOfWindow=data.size()-1;
    int endOfWindow=centerOfWindow;
    int startOfWindow=centerOfWindow;
    double dataCenter=data[data.size()-1];

    //use sliding window to count points within bandwidth distance
    int maxDensity=0;
    int minDensity=-1;
    int steps=0;

    //compute maxDensity in windows of size bandWidth
    while(centerOfWindow>-1)
    {
        if(steps%2==0)
        {
            dataCenter=data[centerOfWindow];
        }
        else
        {
            dataCenter=(data[centerOfWindow+1]+data[centerOfWindow])/2;
        }

        //compute density in window of size bandWidth centered at centerOfWindow
        while(data[endOfWindow]-dataCenter>bandWidth)
        {
            endOfWindow--;
        }
        //now window ends at endOfWindow (inclusive)
        while(startOfWindow>0 && dataCenter-data[startOfWindow-1]<=bandWidth)
        {
            startOfWindow--;
        }
        //now window starts at startOfWindow (inclusive)
        int currDensity=endOfWindow-startOfWindow+1;
        maxDensity=std::max(maxDensity, currDensity); //endOfWindow-startOfWindow+1 is number of points currently within bandWidth
        if(minDensity==-1){minDensity=currDensity;}
        minDensity=std::min(minDensity,currDensity);
        //update of iteration parameters
        if(steps%2==0)
        {
            centerOfWindow--;
        }
        steps++;
    }

   
    //reset iteration parameters
    centerOfWindow=data.size()-1;
    endOfWindow=centerOfWindow;
    startOfWindow=centerOfWindow;
    dataCenter=data[data.size()-1];
    steps=0;

    int lastClusterSplit=data.size()-1;
    //now scan again for places in which density is below significance*maxDensity
    //keep track if point with high density in between
    bool cluster=false;
    int count=0;
    while(centerOfWindow>-1)
    {
        if(steps%2==0)
        {
            dataCenter=data[centerOfWindow];
        }
        else
        {
            dataCenter=(data[centerOfWindow+1]+data[centerOfWindow])/2;
        }

        //compute density in window of size bandWidth centered at centerOfWindow
        while(data[endOfWindow]-dataCenter>bandWidth)
        {
            endOfWindow--;
        }
        //now window ends at endOfWindow (inclusive)
        while(startOfWindow>0 && dataCenter-data[startOfWindow-1]<=bandWidth)
        {
            startOfWindow--;
        }
        //now window starts at startOfWindow (inclusive)
        int currDensity=endOfWindow-startOfWindow+1;
        if(!cluster && currDensity>=maxDensity/(significance)){cluster=true;}
        if(cluster && significance*significance*currDensity<=maxDensity && lastClusterSplit-centerOfWindow+1>=minClusterSize*data.size())
        {
            //detected cluster       
            for(int j=lastClusterSplit; j>=centerOfWindow; j--)
            {
                 (*clusterLabel)[data[j]]=currLabel;
            }
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
   
    return clusterLabel;
}
