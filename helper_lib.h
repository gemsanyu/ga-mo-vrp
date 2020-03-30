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
};

struct Population{
  std::vector<thrust::device_vector<int>> kromosom;
};

struct GetEuclideanDistance : public thrust::binary_function<double,double,double>{
  const double x, y;
  GetEuclideanDistance(double _x, double _y) : x(_x), y(_y) {}

  __device__
  double operator () (double xd, double yd){
    return sqrt((xd-x)*(xd-x)+(yd-y)*(yd-y));
  }
};

struct CheckMaxDistance : public thrust::binary_function<double,double,double>{
  const double maxDistance;
  CheckMaxDistance(double _maxDistance) : maxDistance(_maxDistance) {}

  __device__
  double operator () (double distance, double distanceToDepot){
    if (distance+distanceToDepot>maxDistance){
      return INF_DISTANCE;
    } else {
      return distance;
    }
  }
};

struct CheckMaxCap : public thrust::binary_function<int,double,double>{
  const double maxCap;
  CheckMaxCap(double _maxCap) : maxCap(_maxCap) {}

  __device__
  double operator () (int orderSize, double distance){
    if (orderSize>maxCap){
      return INF_DISTANCE;
    } else {
      return distance;
    }
  }
};

struct CheckIsUsed : public thrust::binary_function<bool,double,double> {
  __device__
  double operator () (bool isUsed, double distance){
    if (isUsed){
      return INF_DISTANCE;
    } else {
      return distance;
    }
  }
};

double getEuclideanDistance(double x0, double y0, double x1, double y1);
void readConfig(std::string configFileName, Config &config);
void readData(Config &config, Data &data);
void initPopulation(Population &population, Data &data, Config &config);
void initKromosomRandom(thrust::device_vector<int> &kromosom, int nCust);
void initKromosomGreedy(thrust::device_vector<int> &kromosom, int initialIdx,
  Data &data, Config &config);

#endif
