//clustering.h
//single header file for clustering project, used to link graph.cpp and heatFlow.cpp
/*
TODO:
-)Destructors?!
-)getDistances for graph has to be changed->use labels, done?
-) update constructors, there are a bunch of new member variables for Graph
*/


//include const in many places...


#include <unordered_map>
#include <iostream>
#include <vector>
#include <queue>
#include <math.h>
#include <stdexcept>
#include <limits>

struct Node;
struct Node{
    Node(void); //default constructor, initialize neighbors to be empty, no dimension or coordinates specified, multiplicity set to 1
    Node(int); //initializes with given mult_;
    Node(int, int=1);  //initializes with given dimension d_, mult_=1;
    Node(const std::vector<double>&, int=1);   //initializes with given neighbors and multiplicity
    Node(const std::vector<double>&,const std::unordered_map<Node*,double>&, int=1);   //initializes with given coordinates,neighbors, and multiplicity
    Node(const std::unordered_map<Node*,double>&, int=1);  //initialized with given coordinates and multiplicity
    std::unordered_map<Node*, double > neighbors_;  //map of neighboring nodes
    int mult_; //multiplicity of node
    int d_; //dimension of data point
    std::vector<double> coordinates_;   //vector of coordinates
};

class Graph{
public:
    Graph(void); //default constructor, initializes emptyy graph
    Graph(std::vector<Node*>); //initializes graph with given vertices
    Graph(Graph&); //copy constructor
    Graph& operator=(Graph&); //copy assignment
    ~Graph(void); // destructor, let's check later if we need to write a custom one

    //basic functions accessing private member variables

    int size (void) const; //returns number of vertices
    Node* getVertex(int) const; //accesses ith vertex if i is in range;
    std::vector<Node*>* getVertices(void); //returns vector of vertices
    double getDistance(Node*,Node*);    //returns distance between two vertices
    int getLabel(Node*);    //returns i for this->vertices_[i]
    bool isConnected(void); //returns true if graph is connected
    int connectedCompCount(void); //returns number of connected components
    int connectedCompLabel(Node*);

    //basic functions changing private member variables

    void insert(Node*);    //insert a new node and edges that come with it
    void insert(Node*,Node*, double); //insert a new, undirected edge with given weight
    void remove(Node*,Node*);   //removes an edge
    void computeDistances(void);    //computes all distances in graph if not done previously


    //functions related to heat diffusion on graph

    //computes laplacian of a function, bool checks if function is well-defined
    std::unordered_map<Node*, double>laplacian(const std::unordered_map<Node*, double>&, bool&); 
    //computes on step in Euler iteration for heat equation
    void heatIterationStep(double timeScale, std::unordered_map<Node*,double>&, bool&);
    //computes as many iterations of heatIterationStep as necessary to get up to time
    void heatDiffusion(double time, double timeScale, std::unordered_map<Node*,double>& func, bool& admissible);
   
        //std::unordered_map<std::pair<Node*,Node*>, double> getDistances(void); //this won't work because pair is not hashable by default
private:
    std::vector<Node*> vertices_;   //vector of vertices comprising graph

    bool labeled_;  //true if labels assigned
    std::unordered_map<Node*, int> label_; //label_[vertices_[i]]=i
    void labelNodes(void); //labels nodes if not done previously


    bool distanceMatrixInitialized_;    //true if distances_ initialized
    bool distancesComputed_;            //true if distances_ filled with actual distances
    std::vector<std::vector<double>> distances_;    //distances_[i][j]=dist(vertices_[i],vertices_[j])

    bool connected_;    //true if graph is connected
    bool connectedChecked_;  //true if we have checked for connectedness
    int connectedCompCount_; //number of connected components
    std::unordered_map<Node*,int> componentLabel_; //labels connected components

    void checkConnected_(void); //checks for connectedness
    void helperDFS_(Node*, int, std::unordered_map<Node*,bool>&, int&);
    void labelComponents_(void);
   
    //std::unordered_map<std::vector<Node*>, double> function_;  
        
};

//kNN functions

Graph* getPrekNN(const std::vector<std::vector<double>>&, int); //returns directed kNN
Graph* getPrekNN(Graph*, int); //returns directed kNN

Graph* getMaxkNNHelper(Graph*); 
Graph* getMaxkNN(const std::vector<std::vector<double>>&, int); //returns max kNN 
Graph* getMaxkNN(Graph*, int); //returns max kNN 

Graph* getMinkNNHelper(Graph);
Graph* getMinkNN(Graph*, int); //returns min kNN 
Graph* getMinkNN(const std::vector<std::vector<double>>&, int); //returns min kNN 


//clustering functions

Node* selectStartNode(const Graph*,const std::unordered_map<Node*,int>*);
std::vector<int> heatClustering(const std::vector<std::vector<double>>&, double minClusterSize=0.1); //returns labels of cluster on data
std::unordered_map<Node*,int> heatClustering(Graph*, double minClusterSize=0.1, double concentrationRadius=1, bool=false);   //returns labels of cluster on data
std::unordered_map<double,int>* oneDimensionalClustering(std::vector<double>&, double, double, double);




