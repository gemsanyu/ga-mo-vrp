#include<algorithm>
#include<iomanip>
#include<iostream>
#include<omp.h>
#include<time.h>

#include"ga_lib.h"
#include"helper_lib.h"
#include"crossover_mutation.h"
using namespace std;

int main(int argc, char **argv){
  srand(time(NULL));
  // ios_base::sync_with_stdio(false);

  char *configFileName = argv[1];
  Config *config = readConfig(configFileName);
  config->fileName = argv[2];
  OrderData *orderData = new OrderData;
  readOrderData(config, orderData);

  // cout << "Using Data : " << config->fileName  <<"\n";
  // cout << "N Customer : " << config->nCust  <<"\n";
  // cout << "Max Distance : " << config->maxDist  <<"\n";
  // cout << "Max Cap : " << config->maxCap  <<"\n";
  // cout << fixed << setprecision(5)<< orderData->depot.x<<" "<<orderData->depot.y<<"\n";
  // cout << "Max Iteration : " << config->maxIter  <<"\n";
  // cout << "convergence threshold : " << config->threshold  <<"\n";
  // cout << "N : " << config->N  <<"\n";
  // cout << "NP : " << config->NP<<"\n";
  // cout << "pc : " << config->pc<<"\n";
  // cout << "pm : " << config->pm<<"\n";
  //
  // cout<<"init_time,iteration,total_time,time_per_iteration,fitness,route_count,total_distance\n";

  double start,end;
  start = omp_get_wtime();
  /*
    Initializing Population of N individu
    initialize by random shuffes and greedy
    but the greedy is not implemented yet :D
    after initialization, evaluate and then sort by fitnes value
  */
  vector<Individu*> population;
  for(int i=0;i<config->N;i++){
    int* kromosom = create1DArrayInt(config->nCust);
    // int *kromosom;
    // cudaMallocHost(&kromosom, config->nCust*sizeof(int));
    if (i%2==0){
      initIndividuRandom(config->nCust, kromosom);
    } else {
      initIndividuGreedy(config, orderData, kromosom);
    }
    Individu* newIdv = new Individu;
    newIdv->kromosom = kromosom;
    decodeKromosom(config, newIdv->kromosom, orderData, &newIdv->routeSet);
    calculateFitness(newIdv);
    population.push_back(newIdv);
  }
  sort(population.begin(), population.end(), cmpIndividuFitness);
  end = omp_get_wtime();
  double initTime = end-start;

  /*
    Prepare the memory in cuda device
    create a minstd_rand object to act as our source of randomness
    create a uniform_int_distribution to produce ints from [0, nCust)
  */
  int *d_nCust;
  cudaMalloc(&d_nCust, sizeof(int));
  cudaMemcpy(d_nCust, &config->nCust, sizeof(int), cudaMemcpyHostToDevice);

  int maxOffSize = config->NP * (config->NP-1);
  vector<cudaStream_t> cudaStreams(maxOffSize);
  vector<int*> odAs(maxOffSize), odBs(maxOffSize);
  vector<int*> mutAs(maxOffSize), mutBs(maxOffSize);
  vector<int*> kromosomAs(maxOffSize), kromosomBs(maxOffSize);
  vector<int*> kromosomOffs(maxOffSize);

  for(int i=0;i<maxOffSize;i++){
    cudaStreamCreate(&cudaStreams[i]);
    cudaMalloc(&odAs[i], sizeof(int));
    cudaMalloc(&odBs[i], sizeof(int));
    cudaMalloc(&mutAs[i], sizeof(int));
    cudaMalloc(&mutBs[i], sizeof(int));
    cudaMalloc(&kromosomAs[i], config->nCust*sizeof(int));
    cudaMalloc(&kromosomBs[i], config->nCust*sizeof(int));
    cudaMalloc(&kromosomOffs[i], config->nCust*sizeof(int));
  }
  int *threadCounts = create1DArrayInt(maxOffSize);
  bool *isMuts = create1DArrayBool(maxOffSize);

  /*
    Start the GA
    for MaxIter
    or until
    the difference between current bestFitness and last bestFitness
    is less than the threshold
    for 100 consecutive iterations
  */
  int sameFitnessCount = 0;
  Individu bestIndividu = *population[0];
  double bestFitness = bestIndividu.fitnessValue;
  // cout<<0<<" "<<bestFitness<<" "<<bestIndividu.routeCount<<" "<<bestIndividu.totalDist<<"\n";
  int iter;
  for (iter=0;iter<config->maxIter && sameFitnessCount<500;iter++){
    /*
      spinning roulette wheel
      then delete the infertile parent based on PC
    */
    vector<int> rwResult = spinRouletteWheel(&population, config->NP);
    vector<int> parentsIdx;
    for(int p=0;p<rwResult.size();p++){
      double r = (double) rand()/RAND_MAX;
      if (r>config->pc){
        continue;
      }
      parentsIdx.push_back(rwResult[p]);
    }

    /*
      generating offsprings
      with complete pairs of the chosen parents doing odx crossover
    */
    int ofIdx = 0;
    int chosenParentSize=parentsIdx.size();
    int offSize = chosenParentSize * (chosenParentSize -1);
    for(int p1=0;p1<chosenParentSize;p1++){
      int pIdx1 = parentsIdx[p1];
      int *kromosomP1 = population[pIdx1]->kromosom;
      for(int p2=p1+1;p2<chosenParentSize;p2++, ofIdx+=2){
        int pIdx2 = parentsIdx[p2];
        int *kromosomP2 = population[pIdx2]->kromosom;
        cudaMemcpyAsync(kromosomAs[ofIdx], kromosomP1,config->nCust*sizeof(int),
          cudaMemcpyHostToDevice, cudaStreams[ofIdx]);
        cudaMemcpyAsync(kromosomBs[ofIdx], kromosomP2, config->nCust*sizeof(int),
          cudaMemcpyHostToDevice, cudaStreams[ofIdx]);
        cudaMemcpyAsync(kromosomAs[ofIdx+1], kromosomBs[ofIdx],config->nCust*sizeof(int),
          cudaMemcpyDeviceToDevice, cudaStreams[ofIdx+1]);
        cudaMemcpyAsync(kromosomBs[ofIdx+1], kromosomAs[ofIdx], config->nCust*sizeof(int),
          cudaMemcpyDeviceToDevice, cudaStreams[ofIdx+1]);
      }
    }

    for(int ofIdx=0;ofIdx<offSize;ofIdx++){
      //order-crossover and mutation points
      int a=rand()%config->nCust;
      int b=rand()%config->nCust;
      if (a>b){
        int c=b;
        b=a;
        a=c;
      }
      cudaMemcpyAsync(odAs[ofIdx], &a, sizeof(int), cudaMemcpyHostToDevice,
        cudaStreams[ofIdx]);
      cudaMemcpyAsync(odBs[ofIdx], &b, sizeof(int), cudaMemcpyHostToDevice,
        cudaStreams[ofIdx]);

      double r = rand()/(double)RAND_MAX;
      if (r<config->pm){
        isMuts[ofIdx]=true;
        a=rand()%config->nCust;
        b=rand()%config->nCust;
        if (a>b){
          int c=b;
          b=a;
          a=c;
        }
        cudaMemcpyAsync(mutAs[ofIdx], &a, sizeof(int), cudaMemcpyHostToDevice,
          cudaStreams[ofIdx]);
        cudaMemcpyAsync(mutBs[ofIdx], &b, sizeof(int), cudaMemcpyHostToDevice,
          cudaStreams[ofIdx]);
        threadCounts[ofIdx]=(b-a+1)/2;
      } else {
        isMuts[ofIdx]=false;
      }
    }

    for(int ofIdx=0;ofIdx<offSize;ofIdx++){
      size_t sharedMemory = config->nCust*sizeof(bool);
      orderCrossover<<< 1, config->nCust, sharedMemory, cudaStreams[ofIdx] >>>(
        d_nCust,
        kromosomAs[ofIdx],
        kromosomBs[ofIdx],
        kromosomOffs[ofIdx],
        odAs[ofIdx],
        odBs[ofIdx]
      );
    }

    for(int ofIdx=0;ofIdx<offSize;ofIdx++){
      if (isMuts[ofIdx]){
        rsMutation<<< 1, threadCounts[ofIdx], 0, cudaStreams[ofIdx] >>>(
          kromosomOffs[ofIdx],
          mutAs[ofIdx],
          mutBs[ofIdx]
        );
      }
    }

    cudaDeviceSynchronize();
    // order crossover & mutation
    for(int i=0;i<offSize; i++){
      int *kromosomOff = create1DArrayInt(config->nCust);
      // int *kromosomOff;
      // cudaMallocHost(&kromosomOff, config->nCust*sizeof(int));
      cudaMemcpy(kromosomOff, kromosomOffs[i],
        config->nCust*sizeof(int), cudaMemcpyDeviceToHost);
      Individu* offspring = new Individu();
      offspring->kromosom = kromosomOff;
      decodeKromosom(config, offspring->kromosom, orderData, &offspring->routeSet);
      calculateFitness(offspring);
      population.push_back(offspring);
    }

    /*
      selection NSGA2
    */
    vector<Individu*> newPopulation = selectionNSGA2(config,&population);
    population = newPopulation;
    vector<Individu*>().swap(newPopulation);


    /*
      note best fitness
      and check for convergence
    */
    sort(population.begin(), population.end(), cmpIndividuFitness);
    double currentBestFitness = population[0]->fitnessValue;
    if (abs(bestFitness-currentBestFitness)<config->threshold){
      sameFitnessCount++;
    } else {
      sameFitnessCount=0;
    }

    if (currentBestFitness>bestFitness){
      bestFitness = currentBestFitness;
      bestIndividu = *population[0];
    }
    // cout<<iter+1<<" "<<bestFitness<<" "<<bestIndividu.routeCount<<" "<<bestIndividu.totalDist<<"\n";
  }
  end = omp_get_wtime();
  double totalTime = end-start;

  // cout<<"Obtained PF\n";
  // for(int i=0;i<population.size();i++){
  //   cout<<fixed<<setprecision(8)<<"solution-"<<i<<": "<<population[i]->routeCount<<" "<<population[i]->totalDist<<"\n";
  // }
  // cout<<fixed<<setprecision(8)<<"Initial Population Generation Time: "<<initTime<<"\n";
  // cout<<"GA stopped at iteration-"<<iter<<"\n";
  // cout<<fixed<<setprecision(8)<<"Total Time: "<<totalTime<<"\n";
  // cout<<fixed<<setprecision(8)<<"Mean Time per iteration: "<<totalTime/(double)iter<<"\n";
  cout<<fixed<<setprecision(8)<<initTime<<",";
  cout<<iter<<",";
  cout<<fixed<<setprecision(8)<<totalTime<<",";
  cout<<fixed<<setprecision(8)<<totalTime/(double)iter<<",";
  cout<<fixed<<setprecision(8)<<bestFitness<<","<<bestIndividu.routeCount<<","<<bestIndividu.totalDist<<"\n";

  /*
    freeing memory
  */
  for(int i=0;i<population.size();i++){
    delete[] population[i]->kromosom;
    delete population[i];
  }
  vector<Individu*>().swap(population);

}
