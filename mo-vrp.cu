#include <iomanip>
#include <iostream>
#include <omp.h>

#include "helper_lib.h"

int main(int argc, char **argv){
  char *configFileName = argv[1];
  Config config;
  readConfig(configFileName, config);
  config.fileName = argv[2];
  Data data;
  readData(config, data);

  Population population;
  thrust::device_vector<int> routeCount(config.N);
  thrust::device_vector<double> totalDist(config.N);
  thrust::device_vector<double> fitnessValue(config.N);

  /*
    init population
  */
  for(int i=0;i<config.N;i++){
    thrust::device_vector<int> newKromosom(config.nCust);
    population.kromosom.push_back(newKromosom);
  }
  initPopulation(population, data, config);
  // decodeKromosom(&population, &data, &config, routeCount, totalDist, fitnessValue);

  // for(int i=0;i<config.N;i++){
  //   for(int j=0;j<config.nCust;j++){
  //     std::cout<<population.kromosom[i][j]<<" ";
  //   }
  //   std::cout<<"\n";
  // }


}
