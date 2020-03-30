#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <random>
#include <vector>


#include"helper_lib.h"

double getEuclideanDistance(double x0, double y0, double x1, double y1){
  return sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
}

void readConfig(std::string configFileName, Config &config){
  std::ifstream configFile(configFileName);

  configFile >> config.maxIter >> config.threshold >>  config.N;
  configFile >> config.NP >> config.pc >> config.pm;
  configFile.close();
}

void readData(Config &config, Data &data){
  std::ifstream dataFile(config.fileName);

  dataFile >> config.nCust >> config.maxDist >> config.maxCap;
  dataFile >> data.depotX >> data.depotY;
  for(int i=0;i<config.nCust;i++){
    double x, y;
    int orderSize;
    dataFile >> x >> y >> orderSize;
    data.customers.x.push_back(x);
    data.customers.y.push_back(y);
    data.customers.orderSize.push_back(orderSize);
  }
  dataFile.close();
}

void initPopulation(Population &population, Data &data, Config &config){
  thrust::default_random_engine randEng;
  thrust::uniform_real_distribution<double> uniDist(0,1);
  randEng.discard(time(NULL));

  for(int i=0;i<config.N;i++){
    if(i%2==0){
      initKromosomRandom(population.kromosom[i], config.nCust);
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

  thrust::device_vector<double> distancesToNextCust(config.nCust);
  thrust::device_vector<double> distancesCustToDepot(config.nCust);
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

  /*
    starting the greedy routing
    first prepare the distance from every customer to depot
  */
  thrust::transform(
    data.customers.x.begin(),
    data.customers.x.end(),
    data.customers.y.begin(),
    distancesCustToDepot.begin(),
    GetEuclideanDistance(depotX, depotY)
  );

  while(servedCustomerCount<config.nCust){
    /*
      finding the nearest customer
      1. calculate distance from last coord to all customers
      2. check isUsed
      3. check MaxDist
      4. check MaxCap
    */
    thrust::transform(
      data.customers.x.begin(),
      data.customers.x.end(),
      data.customers.y.begin(),
      distancesToNextCust.begin(),
      GetEuclideanDistance(lastX, lastY)
    );

    thrust::transform(
      isUsed.begin(),
      isUsed.end(),
      distancesToNextCust.begin(),
      distancesToNextCust.begin(),
      CheckIsUsed()
    );

    thrust::transform(
      distancesToNextCust.begin(),
      distancesToNextCust.end(),
      distancesCustToDepot.begin(),
      distancesToNextCust.begin(),
      CheckMaxDistance(config.maxDist)
    );

    thrust::transform(
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
    int chosenCustIdx=iter-distancesToNextCust.begin();
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

void initKromosomRandom(thrust::device_vector<int> &kromosom, int nCust){
  thrust::sequence(kromosom.begin(), kromosom.end(), 0, 1);
  thrust::host_vector<double> h_randKey(nCust);
  thrust::generate(h_randKey.begin(), h_randKey.end(), rand);
  thrust::device_vector<double> randKey = h_randKey;
  thrust::sort_by_key(randKey.begin(), randKey.end(), kromosom.begin());
}
