 //main.cpp
#include "clustering.h"


#include<iostream>
#include <iomanip>
<<<<<<< HEAD
#include <random>


int main(){ 

    std::vector<std::vector<double>> data;

    //set pseudo random number generator
    std::random_device rd{};
    std::mt19937 generator{rd()};

    //boxes with random noise
    /*
    //generate 500 data points normally distributed around (0,0), correlation matrix is diag(1,2)
    std::uniform_real_distribution<> uniform_x1(-1.0, 1.0);
    std::uniform_real_distribution<> uniform_y1(-1.0, 1.0);
    std::normal_distribution<double> normal_x1(0.0,0.1);
    std::normal_distribution<double> normal_y1(0.0,0.2);
    int size1=1000;

    for(int i=0; i<size1; i++){
        double x = uniform_x1(generator)+normal_x1(generator);
        double y = uniform_y1(generator)+normal_y1(generator);
        data.push_back({x,y});
    }

    //generate 200 data points normally distributed around (0,0), correlation matrix is diag(1,2)
    std::uniform_real_distribution<> uniform_x2(1.3, 3.0);
    std::uniform_real_distribution<> uniform_y2(0.0, 2.0);
    std::normal_distribution<double> normal_x2(0.0,0.2);
    std::normal_distribution<double> normal_y2(0.0,0.1);
    int size2=500;
    
    for(int i=0; i<size2; i++){
        double x = uniform_x2(generator)+normal_x2(generator);
        double y = uniform_y2(generator)+normal_y2(generator);
        data.push_back({x,y});
    }
    */
    

    //circles with random noise
    std::uniform_real_distribution<double> uniform_interval(0, 1.0);
    std::normal_distribution<double> noise(0,0.1);
    const double PI = atan(1)*4;

    //annulus, inner radius 1, outer radius 2, plus normal noise
    for(int i=0; i<500; i++){
        double angle = uniform_interval(generator);
        double radius = uniform_interval(generator);
        double noise_x = noise(generator);
        double noise_y = noise(generator);
        double x = radius*cos(2*PI*angle) + noise_x;
        double y = radius*sin(2*PI*angle) +  noise_y;
        data.push_back({x,y});
    }

    //annulus, inner radius 2.55, outer radius 2.95, plus normal noise
    for(int i=0; i<1000; i++){
        double angle = uniform_interval(generator);
        double radius = 1.5+uniform_interval(generator)/2;
        double noise_x = noise(generator);
        double noise_y = noise(generator);
        double x = radius*cos(2*PI*angle) + noise_x;
        double y = radius*sin(2*PI*angle) +  noise_y;
        data.push_back({x,y});
    }

   


    

=======


int main(){ 
    
    std::vector<std::vector<double>> data;

>>>>>>> parent of f35a29a (Delete main.cpp)
    //uniform square
    /*
   

    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            for(int k=0; k<10; k++){
                data.push_back({0.5*i+0.1,0.5*j,0.5*k+0.2});
            }
        }
    }

    */
    


    
    

   //bunch of squares
<<<<<<< HEAD
    /*
    for(int i=0; i<40; i++){
        for(int j=0; j<40; j++){
            double x=-i;
            x=0.5*x;
            x-=0.2;
            double y=-j;
            y=0.5*y;
            y-=0.2;
=======
    
    for(int i=0; i<20; i++){
        for(int j=0; j<20; j++){
            double x=-i;
            x=0.5*x;
            x-=0.1;
            double y=-j;
            y=0.5*y;
            y-=0.1;
>>>>>>> parent of f35a29a (Delete main.cpp)
            data.push_back({x,y});  
        }
    }
    
    



   
   
    for(int i=0; i<40; i++){
        for(int j=0; j<40; j++){
            double x=0.25*i;
            double y=0.25*j;
            data.push_back({x,y});   
        }
    }

    
<<<<<<< HEAD
     for(int i=0; i<50; i++){
=======
     for(int i=0; i<100; i++){
>>>>>>> parent of f35a29a (Delete main.cpp)
        for(int j=0; j<50; j++){
            double x=i;
            x=0.1*x;
            x+=15;
            double y=j;
            y=0.1*y;
            y+=15;
            data.push_back({x,y});   
        }
    }
    
    
    

    for(int i=1; i<46; i++){
        double x=9.5+0.12*i;
        double y=9.5+0.12*i;
        data.push_back({x,y});
    }

    
<<<<<<< HEAD
    */
=======
    
>>>>>>> parent of f35a29a (Delete main.cpp)

    //circles
   /*
    int n=50;
    double pi = atan(1)*4;
    for(int i=0; i<n; i++){
        double angle= 2*pi/n*i;
        double x=0.2*cos(angle+0.4);
        double y=0.2*sin(angle+0.4);
        data.push_back({x,y});
    }
    for(int i=0; i<n; i++){
        double angle= 2*pi/n*i;
        double x=0.201*cos(angle+0.1);
        double y=0.201*sin(angle+0.1);
        data.push_back({x,y});
    }
    for(int i=0; i<n; i++){
        double angle= 2*pi/n*i;
        double x=0.201*cos(angle+0.2);
        double y=0.199*sin(angle+0.2);
        data.push_back({x,y});
    }
      for(int i=0; i<n; i++){
        double angle= 2*pi/n*i;
        double x=0.199*cos(angle+0.3);
        double y=0.199*sin(angle+0.3);
        data.push_back({x,y});
    }

*/

    /*
    for(int i=0; i<n; i++){
        double angle= 2*pi/n*i;
        double x=cos(angle);
        x*=2;
        double y=sin(angle);
        y*=2;
        data.push_back({x,y});
    }
    */

    
    
    /*

    for(int i=0; i<10; i++){
        double y=0.2;
        y+=i*0.2;
        data.push_back({0,y});
    }
    */
    

    
    
    

    
    
    //thin line connecting them
    

    /*
    for(int i=1; i<1000; i++){
        double y=-1;
        y-=i*0.001;
        data.push_back({0,y});
    }
    */
/*
    for(int i=0; i<70; i++){
        for(int j=0; j<70; j++){
            double x=i;
            x=0.5*x;
            double y=j;
            y=0.5*y;
            data.push_back({x,y});   
        }
    }
*/
    
   

    

    //data generated


    /*
    Graph* testGraph = getMinkNN(data, 4);

    Node*  source = testGraph->getVertex(0);
    for(const auto& [nb, weight] : source->neighbors_){
        std::cout << weight << std::endl;
    }
    
    if(testGraph->isConnected()){
        std::cout << "testGraph is connected" << "\n"; 
    }
    */

    //Graph* G = getMaxkNN(data,3);


    /*
    for (int i = 0; i <G->size(); i++)
    {
        Node* currNode = G->getVertex(i);
        auto it = currNode->neighbors_.begin();
        std::cout << currNode->coordinates_[0] << "," << currNode->coordinates_[1] << std::endl;
        std::cout << "has neighbors:" << std::endl;
        while(it!=currNode->neighbors_.end()){
            Node * nb = it->first;
            double dist=it->second; 
            std::cout << nb->coordinates_[0] << "," << nb->coordinates_[1] << std::endl;
            std::cout << dist << std::endl;
            it++;
        }
       
    }
    */

    
    
    
    //if(G->isConnected()){std::cout << "initally connected" << std::endl;}
    /*
    std::cout << G->getVertex(0)->coordinates_[0] << "," << G->getVertex(0)->coordinates_[1] << std::endl;
    for(const auto& [key,value]: G->getVertex(0)->neighbors_){
        std::cout << key->coordinates_[0] << "," << key->coordinates_[1] << ":"<< value << std::endl;
    }    
    std::cout <<"number of components:" << G->connectedCompCount() << std::endl;
      exit(42);
    */
    //std::unordered_map<Node*,int> heatLabels=heatClustering(data, false);

    
    std::vector<int> heatLabels = heatClustering(data);
    std::unordered_map<int,int> counts;
    for (int i = 0; i < data.size(); i++)
    {
        if(counts.count(heatLabels[i])==0){
            counts[heatLabels[i]]=0;
        }
        counts[heatLabels[i]]++;
        std::cout << "(" << data[i][0] << "," <<data[i][1] << "," 
       // << data[i][2] 
        << ")"
        << ":" << heatLabels[i] << "\n";
    }

    for(const auto& [key, value] : counts){
        if(value>2){
<<<<<<< HEAD
        std::cout << key << ":" << value << std::endl;
=======
        std::cout <<value << std::endl;
>>>>>>> parent of f35a29a (Delete main.cpp)
        }
    }

    /*
    std::cout << "HeatLabels: \n";
    std::unordered_map<int,int> counts;
    for(const auto& [key,value]: heatLabels){
        std::cout <<"(" 
                << key->coordinates_[0]
                << ","
                << key->coordinates_[1]
                << "): "
                << value 
                << std::endl;
        if(counts.count(value)==0){
            counts[value]=0;
        }
        counts[value]++;
    }
    

    
    for(const auto& [key,value]: counts){
        if(value>1){
        std::cout <<  value << "\n";
        }

    }
    */
    
   


    return 0;
}



//OLD TEST CODE
    /*
    TEST LAPLACIAN
    
    
    
    Node* a = new Node;
    Node* b = new Node;
    Node* c = new Node;
    Node* d = new Node;
    Node* e = new Node;
    Node* f = new Node;

    std::vector<Node*> vertices{a,b,c,d,e,f};
    Graph G(vertices);

    G.insert(a,b,1);
    G.insert(b,c,2);
    G.insert(c,d,3);
    G.insert(d,e,4);
    G.insert(a,d,1);
    G.insert(a,e,2);
    G.insert(a,f,12);



    std::unordered_map<Node*,double> func;
    for (int i = 0; i < 6; i++)
    {
        func[G.getVertex(i)]=0;
    }
    func[a]=1;
    bool admissible=false;
    for(const auto& [key,value] : func){
            std::cout << value <<',';
    }
    std::cout << std::endl;
    for(int i=0; i<10; i++){
        G.heatDiffusion(0.1, 0.001, func, admissible);
        for(const auto& [key,value] : func){
            std::cout << value <<',';
        }
        std::cout << std::endl;
    }
    */

    /*
    std::vector<std::vector<double>> Point;
    for (int i = 0; i < 10; i++)
    {
        Point.push_back({0,1});
    }
    Point.push_back({0,2});

    Graph* PointGraph = getMinkNN(Point,11);
    for (int i = 0; i < PointGraph->size(); i++){
        int temp = PointGraph->getVertex(i)->mult_;
        std::cout << temp << "\n";
    }

    std::unordered_map<Node*, double > Distribution;
    for (int i = 0; i < PointGraph->size(); i++)
    {
        Distribution[PointGraph->getVertex(i)]=0;
    }
    Distribution[PointGraph->getVertex(0)]=11;
    bool admissible=false;
    while(Distribution[PointGraph->getVertex(0)]>2){
        PointGraph->heatDiffusion(0.1,0.001,Distribution,admissible);
        std::cout <<Distribution[PointGraph->getVertex(0)] << ", " << Distribution[PointGraph->getVertex(1)] << "\n";

    }
    */



       
   

   /*


    int k=4;
    Graph* prekNN=getPrekNN(data,k);
    Graph* minkNN=getMinkNN(data,k);

    double check=10000;
    for (int i = 0; i < minkNN->size(); i++)
    {
        Node* currNode= minkNN->getVertex(i);
        for(const auto& [key,value] : currNode->neighbors_){
            check=std::min(check,value);
        }
    }
    std::cout << check;
    
    //Node* source = selectStartNode(minkNN);
    //std::cout << source->coordinates_[0] << "," << source->coordinates_[1] << std::endl;
    std::unordered_map<Node*,int> heatLabels=heatClustering(minkNN);
    std::cout<<heatLabels.size();
    
   
    //
    
    for(const auto& [key,value] : heatLabels){
        std::cout << key->coordinates_[0] << "," << key->coordinates_[1] << "\n";
    }
    
    for(int i=0; i<10;i++){
        for(int j=0; j<10; j++){
            int ind = i*10+j;
            std::cout <<prekNN->getVertex(ind)->neighbors_.size() << ",";
        }
        std::cout << std::endl;
    }
    
    std::cout << std::endl;
    for(int i=0; i<10;i++){
        for(int j=0; j<10; j++){
            int ind = i*10+j;
            std::cout <<minkNN->getVertex(ind)->neighbors_.size() << ",";
        }
        std::cout << std::endl;
    }
    

    if(minkNN->isConnected()){
        std::cout << "kNN is connected";
    }
    else{
        std::cout << "kNN is not connected" << "\n";
        std::cout << "kNN has " << minkNN->connectedCompCount() << " components";
    }
    
    

    
    Graph* KNN = getMinkNN(data, k);
    int comp=10000;
    for(int i=0; i<KNN->size();i++){
        int temp=KNN->getVertex(i)->neighbors_.size();
        comp=std::min(comp, temp);
        if(comp==0){
            double x=KNN->getVertex(i)->coordinates_[0];
            double y=KNN->getVertex(i)->coordinates_[1];
            std::cout << x << "," << y << std::endl;
            break;
        }
    }

    std::cout << comp << std::endl;
    std::cout << KNN->connectedCompCount() << std::endl;

    

    

    Graph* kNNFromData= getMinkNN(data, 358);
    if( kNNFromData->isConnected()){
        std::cout << "kNN is connected" << std::endl;;
    }
    else{
        std::cout << "kNN is not connected" << std::endl;
        std::cout << "kNN has " << kNNFromData->connectedCompCount() << " components" << std::endl;
    }
    

    
    for(int i=-9; i<10; i++){
        for(int j=-9; j<10; j++){
            double x=i;
            x/=1000;
            x+=1;
            double y=j;
            y/=1000;
            y+=1;
            data.push_back({x,y});   
        }
    }

    for(int i=-9; i<10; i++){
        for(int j=-9; j<10; j++){
            double x=i;
            x/=1000;
            x+=2;
            double y=j;
            y/=1000;
            y+=2;
            data.push_back({x,y});   
        }
    }
    

   
    
    
    

    
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 50; j++)
        {   
            double x=i*i;
            double y=j;
            double z=i*j;
            data.push_back({x,y,z});
        }   
    }
    
    
    Graph* kNNFromData= getMaxkNN(data, 2);
    //std::vector<Node*>* vertices = kNNFromData->getVertices();
    if( kNNFromData->isConnected()){
        std::cout << "kNN is connected";
    }
    else{
        std::cout << "kNN is not connected";
    }
    
   

    
    
    Node* origin=kNNFromData->getVertex(0);
    Node* anotherNode=kNNFromData->getVertex(1);
    
    for (int i = 0; i < 3; i++)
    {
        std::cout<< origin->coordinates_[i];
    }
    std::cout<< std::endl;

    for( const std::pair<Node*, double>& nbhd : origin->neighbors_){
        Node* nb=nbhd.first; 
        for (int j = 0; j < 3; j++)
        {
            std::cout<< nb->coordinates_[j] << ' ';
        }
        std::cout << std::endl;
    
    }

    for (int i = 0; i < 3; i++)
    {
        std::cout<< anotherNode->coordinates_[i];
    }
    std::cout<< std::endl;

    for( const std::pair<Node*, double>& nbhd : anotherNode->neighbors_){
        Node* nb=nbhd.first; 
        for (int j = 0; j < 3; j++)
        {
            std::cout<< nb->coordinates_[j] << ' ';
        }
        std::cout << std::endl;
    
    }
    
    
    
    
    Node* a = new Node;
    Node* b = new Node;
    Node* c = new Node;
    Node* d = new Node;
    Node* e = new Node;
    Node* f = new Node;

    std::vector<Node*> vertices{a,b,c,d,e,f};
    Graph G(vertices);

    G.insert(a,b,1);
    G.insert(b,c,2);
    G.insert(c,d,3);
    G.insert(d,e,4);
    G.insert(a,d,1);
    G.insert(a,e,2);
    G.insert(a,f,12);

    if(G.isConnected()){
        std::cout << "Graph G is connected" << std::endl;
    }
    else{
        int c = G.connectedCompCount();
        std::cout << c << std::endl;
    }

    


   Graph* H=&G;

   int k=2;

   Graph* kNN=getkNN(H,k);

   std::vector<Node*>* kNNVertices=kNN->getVertices();

   for(int i=0; i<5 ; i++){
        for(auto& [key,value] : (*kNNVertices)[i]->neighbors_){
            std::cout << kNN->getLabel(key) << " ";
        }
        std::cout << std::endl;
   }
    


    */


