#include"ga_lib.h"

__global__ void orderCrossover(int *nCust, int* kromosomA, int* kromosomB, int* kromosomOff,
int *odA, int *odB){

  /*
    copy parentA segment (a,b) to offspring segment(a,b)
  */
  bool *genExistFlag = (bool*) malloc((*nCust)*sizeof(bool));
  for(int i=0;i<(*nCust);i++){
    genExistFlag[i]=false;
  }

  for (int c=*odA;c<=*odB;c++){
    int custID = kromosomA[c];
    kromosomOff[c] = custID;
    genExistFlag[custID] = true;
  }

  /*
    and then add parentB's gens
    not yet contained by the offspring
  */
  int ofIdx=((*odB)+1)%(*nCust);
  for (int genBIdx=((*odB)+1)%(*nCust);ofIdx<(*odA) || ofIdx>(*odB);genBIdx = (genBIdx+1)%(*nCust)){
    int gen = kromosomB[genBIdx];
    if (genExistFlag[gen]){
      continue;
    }
    kromosomOff[ofIdx]=gen;

    genExistFlag[gen]=true;
    ofIdx = (ofIdx+1)%(*nCust);
  }
  free(genExistFlag);
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
