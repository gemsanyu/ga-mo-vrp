#ifndef CROSSOVER_H
#define CROSSOVER_H

using namespace std;

__global__ void orderCrossover(int *nCust, int* kromosomA, int* kromosomB, int* kromosomOff,
  int *odA, int *odB);

__global__ void rsMutationPar(int* kromosomOff, int *mutA, int *mutB);

#endif
