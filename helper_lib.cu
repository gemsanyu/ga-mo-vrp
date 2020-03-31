#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/extrema.h>
#include <thrust/gather.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sort.h>
#include <random>
#include <vector>

#include "helper_lib.h"

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

  /*
    from the custs coord and depot coord
    we can pre-compute distances among customers and cust-depot
  */
  thrust::device_vector<double> distancesToCust(config.nCust*config.nCust);
  thrust::device_vector<double> distancesToDepot(config.nCust);

  thrust::transform(
    data.customers.x.begin(),
    data.customers.x.end(),
    data.customers.y.begin(),
    distancesToDepot.begin(),
    GetEuclideanDistance(data.depotX, data.depotY)
  );

  for(int i=0;i<config.nCust;i++){
    double custX=data.customers.x[i];
    double custY=data.customers.y[i];

    thrust::transform(
      data.customers.x.begin(),
      data.customers.x.end(),
      data.customers.y.begin(),
      distancesToCust.begin()+i*config.nCust,
      GetEuclideanDistance(custX, custY)
    );
  }

  data.distancesToCust = distancesToCust;
  data.distancesToDepot = distancesToDepot;
}

void sortPopulationByFitness(Population &population, Config const &config){

  /*
    prepare space for sorted values
  */
  thrust::device_vector<double> sortedTotalDist(config.nCust);
  thrust::device_vector<int> sortedRouteCount(config.nCust);

  /*
    initialize indices vector to [0,1,2,..]
  */
  thrust::counting_iterator<int> iter(0);
  thrust::device_vector<int> indices(config.N);
  thrust::copy(iter, iter + indices.size(), indices.begin());
  
  /*
    first sort the keys and indices by the keys
  */
  thrust::sort_by_key(
    population.fitnessValue.begin(),
    population.fitnessValue.end(),
    indices.begin(),
    thrust::greater<double>()
  );

  /*
    Now reorder totalDist, routeCount, fitnessValue and kromosom
    using the sorted indices
  */
  thrust::gather(
    indices.begin(),
    indices.end(),
    population.totalDist.begin(),
    sortedTotalDist.begin()
  );
  population.totalDist = sortedTotalDist;

  thrust::gather(
    indices.begin(),
    indices.end(),
    population.routeCount.begin(),
    sortedRouteCount.begin()
  );
  population.routeCount = sortedRouteCount;

  std::vector<thrust::device_vector<int>> sortedKromosom(config.N);
  for(int i=0;i<config.N;i++){
    sortedKromosom[i]=population.kromosom[indices[i]];
  }
  population.kromosom = sortedKromosom;
}
