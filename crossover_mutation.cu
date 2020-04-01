#include "helper_lib.h"

void crossoverMutation(Population population, Population &offspring,
  thrust::device_vector<int> parentsIdx, Config config){

    int *d_nCust;
    cudaMalloc(&d_nCust, sizeof(int));
    cudaMemcpy(d_nCust, &config.nCust, sizeof(int), cudaMemcpyHostToDevice);
    int offCount = offspring.kromosom.size();
    int parentCount = parentsIdx.size();

  /*
    Crossover and mutation
    preparing
    1. crossover point
    2. mutation point
    3. probability of mutation per kromosom
  */
  thrust::device_vector<int*> d_kromosomOffs(offCount);
  for(int i=0;i<offCount;i++){
    d_kromosomOffs[i]=thrust::raw_pointer_cast(offspring.kromosom[i].data());
  }
  int **pd_kromosomOffs = thrust::raw_pointer_cast(d_kromosomOffs.data());

  thrust::host_vector<int> h_odA(offCount), h_odB(offCount);
  thrust::generate(h_odB.begin(), h_odB.end(), RandInt(config.nCust));
  thrust::transform(
    h_odB.begin(),
    h_odB.end(),
    h_odA.begin(),
    RandIntDynamicUpper()
  );
  thrust::device_vector<int> odA=h_odA;
  thrust::device_vector<int> odB=h_odB;
  int *p_odA = thrust::raw_pointer_cast(odA.data());
  int *p_odB = thrust::raw_pointer_cast(odB.data());

  thrust::host_vector<int> h_mutA(offCount), h_mutB(offCount);
  thrust::generate(h_mutB.begin(), h_mutB.end(), RandInt(config.nCust));
  thrust::transform(
    h_mutB.begin(),
    h_mutB.end(),
    h_mutA.begin(),
    RandIntDynamicUpper()
  );
  
  thrust::device_vector<int> mutA=h_mutA;
  thrust::device_vector<int> mutB=h_mutB;
  int *p_mutA = thrust::raw_pointer_cast(mutA.data());
  int *p_mutB = thrust::raw_pointer_cast(mutB.data());

  thrust::host_vector<double> h_mutProb(offCount);
  thrust::generate(h_mutProb.begin(), h_mutProb.end(), Rand01());
  thrust::device_vector<double> mutProb = h_mutProb;

  thrust::device_vector<int*> d_kromosomAs;
  thrust::device_vector<int*> d_kromosomBs;
  int ofIdx=0;
  for(int i=0;i<parentCount;i++){
    int pIdx1=parentsIdx[i];
    for(int j=0;j<parentCount;j++, ofIdx++){
      if(i==j){
        continue;
      }
      int pIdx2=parentsIdx[j];
      d_kromosomAs.push_back(thrust::raw_pointer_cast(population.kromosom[pIdx1].data()));
      d_kromosomBs.push_back(thrust::raw_pointer_cast(population.kromosom[pIdx2].data()));
    }
  }
  int **pd_kromosomAs = thrust::raw_pointer_cast(d_kromosomAs.data());
  int **pd_kromosomBs = thrust::raw_pointer_cast(d_kromosomBs.data());

  /*
    crossover
  */
  size_t sharedMemory = config.nCust*sizeof(bool);
  orderCrossover<<< offCount, config.nCust, sharedMemory >>>(
    d_nCust,
    pd_kromosomAs,
    pd_kromosomBs,
    pd_kromosomOffs,
    p_odA,
    p_odB
  );

  /*
    mutation
  */
  rsMutationPar<<< offCount, config.nCust >>>(
    pd_kromosomOffs,
    p_mutA,
    p_mutB
  );
  cudaFree(d_nCust);
}

__global__ void orderCrossover(int *nCust, int **kromosomAs, int **kromosomBs,
  int **kromosomOffs, int *odAs, int *odBs){

  int bIdx = blockIdx.x;
  int tIdx = threadIdx.x;
  extern __shared__ bool genExistFlag[];
  genExistFlag[tIdx]=false;
  int *kromosomA = kromosomAs[bIdx];
  int *kromosomB = kromosomBs[bIdx];
  int *kromosomOff = kromosomOffs[bIdx];
  int odA = odAs[bIdx];
  int odB = odBs[bIdx];

  /*
    copy parentA segment (a,b) to offspring segment(a,b)
  */
  if(tIdx >=odA && tIdx<=odB ){
    int custID = kromosomA[tIdx];
    kromosomOff[tIdx] = custID;
    genExistFlag[custID] = true;
  }
  __syncthreads();

  /*
    and then add parentB's gens
    not yet contained by the offspring
  */
  if(tIdx==0){
    int ofIdx=(odB+1)%(*nCust);
    for (int genBIdx=(odB+1)%(*nCust); ofIdx<odA || ofIdx>odB;genBIdx = (genBIdx+1)%(*nCust)){
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

__global__ void rsMutationPar(int **kromosomOffs, int *mutAs, int *mutBs){
  int bIdx = blockIdx.x;
  int tIdx = threadIdx.x;

  int *kromosomOff = kromosomOffs[bIdx];
  int mutA = mutAs[bIdx];
  int mutB = mutBs[bIdx];
  if (tIdx+mutA <= mutB && mutB-tIdx>=mutA){
    int indxMutA = mutA + tIdx;
    int indxMutB = mutB - tIdx;

    //Swapping Algorithm
    int custID = kromosomOff[indxMutA];
    kromosomOff[indxMutA] = kromosomOff[indxMutB];
    kromosomOff[indxMutB] = custID;
  }
}
