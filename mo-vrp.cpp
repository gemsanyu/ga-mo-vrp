#include<iostream>

#include"ga_lib.h"
#include"helper_lib.h"
using namespace std;

int main(int argc, char **argv){
  char *configFileName = argv[1];
  Config config = readConfig(configFileName);
  OrderData orderData = readOrderData(config);

  /*
    Initializing Population of N individu
    initialize by random shuffes
  */
  Individu *population = create1DArrayIndividu(config.N);
  for(int i=0;i<config.N;i++){
    population[i]=initIndividuRandom(config.nCust);
    population[i].routeSet = decodeKromosom(config, population[i].kromosom, orderData);
    calculateFitness(&population[i]);
  }

  /*
    Start the GA
    for MaxIter
    or until the difference between
    the difference between current bestFitness and last bestFitness
    is less than the threshold
    for 100 consecutive iterations
  */
  int sameFitnessCount = 0;
  Individu bestIndividu;
  for (int iter=0;iter<config.maxIter;iter++){
    
  }

}
