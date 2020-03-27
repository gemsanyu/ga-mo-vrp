#include"ga_lib.h"

__global__ void orderCrossover(int *nCust, int* kromosomA, int* kromosomB, int* kromosomOff,
int *odA, int *odB){

  /*
    copy parentA segment (a,b) to offspring segment(a,b)
  */
  extern __shared__ bool genExistFlag[];
  int idx=threadIdx.x;
  int stride=blockDim.x;
  for(int i=idx;i<(*nCust);i+=stride){
    genExistFlag[i]=false;
  }

  int custID = kromosomA[(*odA)+idx];
  kromosomOff[(*odA)+idx] = custID;
  genExistFlag[custID] = true;
  // __syncthreads();

  /*
    and then add parentB's gens
    not yet contained by the offspring
  */
  if(idx==0){
    int ofIdx=((*odB)+1)%(*nCust);
    for (int genBIdx=((*odB)+1)%(*nCust);ofIdx<(*odA) || ofIdx>(*odB);genBIdx = (genBIdx+1)%(*nCust)){
      int gen = kromosomB[genBIdx];
      if (genExistFlag[gen]){
        continue;
      }
      kromosomOff[ofIdx]=gen;
      ofIdx = (ofIdx+1)%(*nCust);
    }
  }
}

__global__ void rsMutationPar(int* kromosomOff, int *mutA, int *mutB){

  int idx=threadIdx.x;
  int indxMutA = *mutA + idx;
  int indxMutB = *mutB - idx;

  //Swapping Algorithm
  int custID = kromosomOff[indxMutA];
  kromosomOff[indxMutA] = kromosomOff[indxMutB];
  kromosomOff[indxMutB] = custID;
}
