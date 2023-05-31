//graph.h


//header file for graph.cpp, node.cpp, knn.cpp, heat_flow.cpp
#ifndef GRAPH_H
#define GRAPH_H


#include <vector>
#include <unordered_map>
#include <iostream>
#include <exception>
#include <memory>
#include <queue>
#include <algorithm>
#include <math.h>
#include <set>


class Node;
//neighbors may not be symmetric

class Node{
public:
//only member access functions
    std::vector<double> getCoordinates(void) const;
private:
    friend class Graph;
    friend class kNNHelper;
    //constructors
    Node(void); //default, initilizes neighbors to be empty, mult_ to 1, dimension to -1, coordinates empty
    Node(const int mult); //initializes with given mult_;
    Node(const int d, int mult=1);  //initializes with given dimension d_, mult_=1;
    Node(const std::vector<double>&, int mult=1);   //initializes with given neighbors and multiplicity
    Node(const std::vector<std::pair<Node*,double>>&, int mult=1);  //initialized with given coordinates and multiplicity
    Node(const std::vector<double>&, const std::vector<std::pair<Node*,double>>&, int mult_=1);   //initializes with given coordinates,neighbors, and multiplicity
    
    //copy,assignment,destructor all defauls
    //beware of what they do!

    //~Node(void);


    /// @brief  adds neighbor of given distance
    /// @param  node to be added as neighbor 
    /// @param  distance to neighbor
    void addNeighbor(Node&, double); //adds neighbor with given weight, if neighbor present updates weight
    void addNeighbor(std::pair<Node*,double>); //calls above function
    void removeNeighbor(Node& node); //removes neighbor

    std::vector<std::pair<Node*,double>>::iterator removeNeighborIt(Node& node);

    bool checkNeighbor(const Node&) const; 


    std::vector<std::pair<Node*,double>> neighbors_; // [neighbor node, dist] 
    //std::unordered_map<Node*,double> neighbors_;  //map of neighbors [neighbor,distance]
    int mult_; //multiplicity of node
    const int d_; //dimension of data point
    std::vector<double> coordinates_;   //vector of coordinates
};

//Graph
//constructors create copies on heap of nodes
//neighbors are symmetric
//at most on edge between nodes
//never adds edges to nodes not contained in graph

class Graph{
public:
    Graph(void); //default constructor, initializes empty graph
    Graph(const std::vector<Node*>& vertices); //initialized with copies of given nodes on stack, this makes dealing with destructor easier

    //copy and assignment
    Graph(const Graph&);
    Graph& operator=(const Graph&);

    //destructor
    ~Graph(void); //deletes all corresponding Nodes

    //member access functions
    int size(void) const;
    Node& getVertex(int) const;
    const std::vector<Node*>& getVertices(void) const;
    int findNode(const Node&) const; //returns index of node if present, otherwise -1
    bool checkNeighbors(const Node&, const Node&, bool nodesOK=false) const;
    int countNeighbors(const Node&) const;
    int countNeighbors(int) const;

    //change member variables
    void insert(const Node&); //insert (copy of) node
    void insert(Node&, Node&, double); //insert edge
    void remove(Node*); //remove and delete node.
    void remove(Node&, Node&);
    //void remove(Node*, Node*); //removes edge between nodes if previously present


    double getMinDistance(void);
    bool checkSym(void);

    /// @brief normalizes distances to newMinDist(=1 by defaul)
    /// @param  newMinDist new minimum distance
    void normalizeDistances(double newMinDist=1);

    bool isConnected(void);
    
    
    /// @brief computes graph laplacian of function
    /// @param vector with ordering same vertices_ 
    std::vector<double> laplacian(const std::vector<double>&);
    void heatIterationStep(std::vector<double>& initialDist, double timeScale, double& low, double& high);
 
    

private:
    std::vector<Node*> vertices_;
    //by constructors those Nodes are always on the heap
    std::vector<std::vector<std::pair<int,double>>> adjacencyLists_; //adjacencyLists[i] contains [j,dist_{ij}] if vertices_[i] and vertices_[j] are nbs
    friend class kNNHelper;
    friend class prekNNHelper;
};

//only required for constructing kNNs


/// @brief computes minimal k nearest neighbor graph of data
/// @param data vector of data points
/// @param k two points are neighbors iff they both among k closest to other
/// @return unq
std::unique_ptr<Graph> getkNN(std::vector<std::vector<double>>& data, int k);
std::unique_ptr<Graph> getkNN(const Graph& G, int k);


#endif //GRAPH_H
