#include "crossover_mutation.h"

__global__ void orderCrossover(int *nCust, int *kromosomA, int *kromosomB,
  int *kromosomOff, int *odA, int *odB){


  int tIdx = threadIdx.x;

  extern __shared__ bool genExistFlag[];
  genExistFlag[tIdx] = false;

  /*
    copy parentA segment (a,b) to offspring segment(a,b)
  */
  if (tIdx >= (*odA) && tIdx <= (*odB)){
    int custID = kromosomA[tIdx];
    genExistFlag[custID]=true;
    kromosomOff[tIdx]=custID;
  }
  
  /*
    and then add parentB's gens
    not yet contained by the offspring
  */
  if (tIdx==0){
    int ofIdx=(*odB+1)%(*nCust);
    for (int genBIdx=(*odB+1)%(*nCust);ofIdx<*odA || ofIdx>*odB;genBIdx = (genBIdx+1)%(*nCust)){
      int gen = kromosomB[genBIdx];
      if (genExistFlag[gen]){
        continue;
      }
      kromosomOff[ofIdx]=gen;
      genExistFlag[gen]=true;
      ofIdx = (ofIdx+1)%(*nCust);
    }
  }

}

__global__ void rsMutation(int *kromosom, int *mutA, int *mutB){

  int tIdx = threadIdx.x;
  /*
    First randomize Mutation-segment points a and b
  */
  int indxMutA = *mutA+tIdx;
  int indxMutB = *mutB-tIdx;

  //Parallel Swapping Algorithm
  int custID = kromosom[indxMutA];
  kromosom[indxMutA] = kromosom[indxMutB];
  kromosom[indxMutB] = custID;
}
