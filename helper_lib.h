#ifndef HLP_LIB_H
#define HLP_LIB_H
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/random.h>
#include <cstring>
#include <vector>
#include <random>

const double INF_DISTANCE=1000000000;

struct Config{
  std::string fileName;
  int nCust;
  int maxDist, maxCap;
  double threshold;
  int maxIter, N, NP;
  double pc, pm;
};

struct Customers{
  thrust::device_vector<double> x,y;
  thrust::device_vector<int> orderSize;
};

struct Data{
  double depotX, depotY;
  Customers customers;
  thrust::device_vector<double> distancesToCust;
  thrust::device_vector<double> distancesToDepot;
};

struct Population{
  std::vector<thrust::device_vector<int>> kromosom;
  thrust::device_vector<int> routeCount;
  thrust::device_vector<double> totalDist;
  thrust::device_vector<double> fitnessValue;
};

struct GetEuclideanDistance : public thrust::binary_function<double,double,double>{
  const double x, y;
  GetEuclideanDistance(double _x, double _y) : x(_x), y(_y) {}

  __device__
  double operator () (double xd, double yd){
    return sqrt((xd-x)*(xd-x)+(yd-y)*(yd-y));
  }
};

struct CalculateFitness : public thrust::binary_function<double,int,double>{
  __device__
  double operator () (double totalDist, int routeCount){
    return 1/(600.0*totalDist + 10000.0*(double)routeCount);
  }
};

struct Rand01{
  __host__
  double operator () (){
    return ((double)rand() / (double)RAND_MAX);
  }
};

struct RandInt{
  int upper;

  RandInt(int _upper): upper(_upper) {}
  __host__
  double operator () (){
    int ret =  rand()%upper;
    if (ret==0){
      ret++;
    }
    return ret;
  }
};

struct RandIntDynamicUpper{
  __host__
  int operator () (int upper){
    return rand()%upper;
  }
};

struct Lesser{
  double r;
  Lesser(double _r) : r(_r) {}

  __device__
  bool operator () (double val){
    return val<=r;
  }
};

struct Greater{
  double r;
  Greater(double _r) : r(_r) {}

  __device__
  bool operator () (double val){
    return val>=r;
  }
};


double getEuclideanDistance(double x0, double y0, double x1, double y1);
void readConfig(std::string configFileName, Config &config);
void readData(Config &config, Data &data);
void initPopulation(Population &population, Data &data, Config &config);
void initKromosomRandom(thrust::device_vector<int> &kromosom, int nCust);
void initKromosomGreedy(thrust::device_vector<int> &kromosom, int initialIdx,
  Data &data, Config &config);
void decodeKromosom(thrust::device_vector<int> kromosom, Data data,
  Config config, int &routeCount, double &totalDist);
void sortPopulationByFitness(Population &population, Config const &config);
void getParentsIdx(Population const &population, Config const &config,
  thrust::device_vector<int> &parentsIdx, int &parentCount);
  void crossoverMutation(Population population, Population &offspring,
    thrust::device_vector<int> parentsIdx, Config config);
__global__ void orderCrossover(int *nCust, int **kromosomAs, int **kromosomBs,
  int **kromosomOffs, int *odAs, int *odBs);
__global__ void rsMutationPar(int **kromosomOffs, int *mutAs, int *mutBs);
#endif
