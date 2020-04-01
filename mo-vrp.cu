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

  /*
    init population
  */
  for(int i=0;i<config.N;i++){
    thrust::device_vector<int> newKromosom(config.nCust);
    population.kromosom.push_back(newKromosom);
  }
  initPopulation(population, data, config);

  /*
    decoding and computing fitness value
  */
  population.totalDist = thrust::device_vector<double>(config.N);
  population.fitnessValue = thrust::device_vector<double>(config.N);
  population.routeCount = thrust::device_vector<double>(config.N);
  for(int i=0;i<config.N;i++){
    int routeCount_t;
    double totalDist_t;
    decodeKromosom(population.kromosom[i], data, config,
    routeCount_t, totalDist_t);
    population.routeCount[i]=routeCount_t;
    population.totalDist[i]=totalDist_t;
  }

  thrust::transform(
    population.totalDist.begin(),
    population.totalDist.end(),
    population.routeCount.begin(),
    population.fitnessValue.begin(),
    CalculateFitness()
  );
  sortPopulationByFitness(population, config);

  /*
    start the GA iteration wohoo !
  */
  double bestFitness = population.fitnessValue[0];
  thrust::device_vector<int> bestKromosom = population.kromosom[0];
  int sameFitnessCount=0;

  for(int t=0;t<config.maxIter && sameFitnessCount<500; t++){
    /*
      parent selection by roulette wheel based on
      fitness value
    */
    thrust::device_vector<int> parentsIdx(config.NP);
    int parentCount;
    getParentsIdx(population, config, parentsIdx, parentCount);

    /*
      crossover and mutation
    */
    Population offspring;
    int offCount=2*parentCount*(parentCount-1);
    for(int i=0;i<offCount;i++){
      thrust::device_vector<int> kromosom(config.nCust,0);
      offspring.kromosom.push_back(kromosom);
    }
    offspring.totalDist = thrust::device_vector<double>(offCount);
    offspring.fitnessValue = thrust::device_vector<double>(offCount);
    offspring.routeCount = thrust::device_vector<double>(offCount);
    crossoverMutation(population, offspring, parentsIdx, config);
  }


  // for(int i=0;i<config.N;i++){
  //   for(int j=0;j<config.nCust;j++){
  //     std::cout<<population.kromosom[i][j]<<" ";
  //   }
  //   std::cout<<"\n";
  // }
  //
  // for(int i=0;i<data.distancesToDepot.size();i++){
  //   std::cout<<data.distancesToDepot[i]<<" ";
  // }
  // std::cout<<"\n";

}
