#include "clustering.h"

#include<iostream>
#include <random>

//sample some data and attempt clustering
void testClustering(){


    //vector of data points
    std::vector<std::vector<double>> data;

    //generate data below

    //set pseudo random number generator
    std::random_device rd{};
    std::mt19937 generator{rd()};

    //probability distributions
    std::uniform_real_distribution<double> uniform_interval(0, 1.0);
    std::normal_distribution<double> noise(0,1);

    const double PI = atan(1)*4;


    //circle with radius 5/6 and some normal noise
    for(int i=0; i<250; i++){
        double angle = uniform_interval(generator);
        double radius = uniform_interval(generator)/1.2;
        double noise_x = noise(generator);
        double noise_y = noise(generator);
        double x = radius*cos(2*PI*angle) + noise_x;
        double y = radius*sin(2*PI*angle) +  noise_y;
        data.push_back({x,y});
    }
    

    //annulus, inner radius 5, outer radius 5.1, plus normal noise
    for(int i=0; i<1000; i++){
        double angle = uniform_interval(generator);
        double radius = 5+uniform_interval(generator)/10;
        double noise_x = noise(generator);
        double noise_y = noise(generator);
        double x = radius*cos(2*PI*angle) + noise_x/2;
        double y = radius*sin(2*PI*angle) +  noise_y;
        data.push_back({x,y});
    }

    //box [7,8]x[0,1] plus noise
    for(int i=0; i<500; i++){
        double x = 7+uniform_interval(generator)+noise(generator)/2;
        double y = uniform_interval(generator)+noise(generator)/15;
        data.push_back({x,y});

    }

    //data generated

    //cluster
    std::vector<int> heatLabels = heatClustering(data);

    //output results
    std::unordered_map<int,int> counts;
    for (int i = 0; i < data.size(); i++)
    {
        if(counts.count(heatLabels[i])==0){
            counts[heatLabels[i]]=0;
        }
        counts[heatLabels[i]]++;
        std::cout << "(" << data[i][0] << "," <<data[i][1]<< ")"
        << ":" << heatLabels[i] << "\n";
    }

    for(const auto& [key, value] : counts){
        if(value>2){
        std::cout << key << ":" << value << std::endl;
        }
    }

    return;

}
