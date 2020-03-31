#ifndef HLP_LIB_H
#define HLP_LIB_H
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/random.h>
#include <cstring>
#include <vector>

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

#endif
