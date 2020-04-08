#ifndef CROSSOVER_MUTATION_H
#define CROSSOVER_MUTATION_H

using namespace std;

__global__ void orderCrossover(int *nCust, int *kromosomA, int *kromosomB,
  int *kromosomOff, int *odA, int *odB);
__global__ void rsMutation(int *kromosomOff, int *mutA, int *mutB);

#endif
