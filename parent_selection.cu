#include <thrust/scan.h>
#include <thrust/find.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/random.h>
#include <thrust/sort.h>
#include <random>
#include <vector>

#include "helper_lib.h"

struct Normalize{
  const double sumFitness;

  Normalize(double _sumFitness) : sumFitness(_sumFitness) {}

  __device__
  double operator () (double fitnessValue){
    return fitnessValue/sumFitness;
  }
};


void getParentsIdx(Population const &population, Config const &config,
  thrust::device_vector<int> &parentsIdx, int &parentCount){

  double sumFitness;
  sumFitness = thrust::reduce(
    population.fitnessValue.begin(),
    population.fitnessValue.end(),
    0.,
    thrust::plus<double>()
  );

  /*
    get probability (normalized fitness values)
    then prefix sum the probability
    to ease searching
  */
  thrust::device_vector<double> probs(config.N); //based on fitness
  thrust::transform(
    population.fitnessValue.begin(),
    population.fitnessValue.end(),
    probs.begin(),
    Normalize(sumFitness)
  );

  thrust::inclusive_scan(
    probs.begin(),
    probs.end(),
    probs.begin()
  );

  /*
    generate random value and spin roulette wheel
  */
  srand(time(NULL));
  thrust::host_vector<double> h_randKey(config.NP);
  thrust::generate(h_randKey.begin(), h_randKey.end(), Rand01());
  thrust::device_vector<double> randKey = h_randKey;
  thrust::device_vector<int> rwParentsIdx(config.NP);

  for(int i=0;i<config.NP;i++){
    thrust::device_vector<double>::iterator iter1, iter2;
    iter1 = probs.begin();
    iter2 = thrust::find_if(
      probs.begin(),
      probs.end(),
      Greater(randKey[i])
    );
    rwParentsIdx[i] = thrust::distance(iter1, iter2);
  }

  /*
    check parents fertility (PC)
    generate random value 1 more time, then copy if fertile
  */
  thrust::device_vector<int> chosenParentsIdx(config.NP);
  thrust::generate(h_randKey.begin(), h_randKey.end(), Rand01());
  randKey = h_randKey;

  thrust::device_vector<int>::iterator finalParentEnd = thrust::copy_if(
    rwParentsIdx.begin(),
    rwParentsIdx.end(),
    randKey.begin(),
    chosenParentsIdx.begin(),
    Lesser(config.pc)
  );

  parentCount = thrust::distance(chosenParentsIdx.begin(), finalParentEnd);
  thrust::copy(
    chosenParentsIdx.begin(),
    chosenParentsIdx.begin()+parentCount,
    parentsIdx.begin()
  );
}
