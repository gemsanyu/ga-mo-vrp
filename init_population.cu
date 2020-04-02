#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/random.h>
#include <thrust/system/cuda/execution_policy.h>
#include <omp.h>
#include <cstring>
#include <vector>

#include "helper_lib.h"


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

struct GenRand{
  thrust::default_random_engine randEng;
  thrust::uniform_real_distribution<double> uniDist;

  GenRand(thrust::default_random_engine _randEng,
    thrust::uniform_real_distribution<double> _uniDist) : randEng(_randEng), uniDist(_uniDist) {}

  __device__
  double operator () (){
    return uniDist(randEng);
  }
};


void initPopulation(Population &population, Data &data, Config &config){
  thrust::default_random_engine randEng;
  thrust::uniform_real_distribution<double> uniDist(0,1);
  randEng.discard(time(NULL));

  #pragma omp parallel for num_threads(config.N)
  for(int i=0;i<config.N;i++){
    if(i%2==0){
      initKromosomRandom(population.kromosom[i], config.nCust, randEng, uniDist);
    } else {
      int initialIdx = uniDist(randEng)*(double)config.nCust;
      initKromosomGreedy(population.kromosom[i], initialIdx, data, config);
    }
  }
}

void initKromosomGreedy(thrust::device_vector<int> &kromosom, int initialIdx,
  Data &data, Config &config){

  double depotX = data.depotX;
  double depotY = data.depotY;

  int maxCap = config.maxCap;
  double maxDist = config.maxDist;
  int servedCustomerCount = 0;

  /*
    starting the greedy routing
    first prepare the distance from every customer to depot
  */
  thrust::device_vector<double> distancesCustToDepot = data.distancesToDepot;
  thrust::device_vector<bool> isUsed(config.nCust, false);

  /*
    get random initial customer
  */
  double lastX = data.customers.x[initialIdx];
  double lastY = data.customers.y[initialIdx];

  double dist = getEuclideanDistance(depotX, depotY, lastX, lastY);
  maxDist = maxDist - dist;
  maxCap = maxCap - data.customers.orderSize[initialIdx];
  kromosom[servedCustomerCount]=initialIdx;
  servedCustomerCount++;
  isUsed[initialIdx]=true;
  int chosenCustIdx = initialIdx;

  while(servedCustomerCount<config.nCust){
    /*
      finding the nearest customer
      1. calculate distance from last coord to all customers
      2. check isUsed
      3. check MaxDist
      4. check MaxCap
    */
    thrust::device_vector<double> distancesToNextCust(
      data.distancesToCust.begin() + chosenCustIdx*config.nCust,
      data.distancesToCust.begin() + (chosenCustIdx+1)*config.nCust
    );

    thrust::transform(
      thrust::device,
      isUsed.begin(),
      isUsed.end(),
      distancesToNextCust.begin(),
      distancesToNextCust.begin(),
      CheckIsUsed()
    );

    thrust::transform(
      thrust::device,
      distancesToNextCust.begin(),
      distancesToNextCust.end(),
      distancesCustToDepot.begin(),
      distancesToNextCust.begin(),
      CheckMaxDistance(config.maxDist)
    );

    thrust::transform(
      thrust::device,
      data.customers.orderSize.begin(),
      data.customers.orderSize.end(),
      distancesToNextCust.begin(),
      distancesToNextCust.begin(),
      CheckMaxCap(config.maxCap)
    );

    thrust::device_vector<double>::iterator iter = thrust::min_element(
      distancesToNextCust.begin(),
      distancesToNextCust.end()
    );
    chosenCustIdx=iter-distancesToNextCust.begin();
    double closestDistance=*iter;

    /*
      check if feasible customer to visit exists
      else go back to depot and start again
    */
    if(closestDistance==INF_DISTANCE){
      maxDist=config.maxDist;
      maxCap=config.maxCap;
      lastX=depotX;
      lastY=depotY;
      continue;
    }

    maxDist = maxDist - closestDistance;
    maxCap = maxCap - data.customers.orderSize[chosenCustIdx];
    lastX = data.customers.x[chosenCustIdx];
    lastY = data.customers.y[chosenCustIdx];
    kromosom[servedCustomerCount]=chosenCustIdx;
    servedCustomerCount++;
    isUsed[chosenCustIdx]=true;
  }
}

void initKromosomRandom(thrust::device_vector<int> &kromosom, int nCust,
  thrust::default_random_engine randEng, thrust::uniform_real_distribution<double> uniDist){
  thrust::sequence(kromosom.begin(), kromosom.end(), 0, 1);
  thrust::device_vector<double> randKey(nCust);
  thrust::generate(randKey.begin(), randKey.end(), GenRand(randEng, uniDist));
  thrust::sort_by_key(randKey.begin(), randKey.end(), kromosom.begin());
}
