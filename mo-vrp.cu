#include<algorithm>
#include<iomanip>
#include<iostream>
#include<omp.h>
#include<time.h>

#include"crossover.h"
#include"ga_lib.h"
#include"helper_lib.h"
using namespace std;

int main(int argc, char **argv){
  srand(time(NULL));
  // ios_base::sync_with_stdio(false);

  char *configFileName = argv[1];
  Config *config = readConfig(configFileName);
  config->fileName = argv[2];
  OrderData *orderData = new OrderData;
  readOrderData(config, orderData);

  int* d_nCust;
  cudaMalloc(&d_nCust, sizeof(int));
  cudaMemcpy(d_nCust, &config->nCust, sizeof(int), cudaMemcpyHostToDevice);

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
  cout<<"init_time,iteration,total_time,time_per_iteration,fitness,route_count,total_distance\n";

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
    vector<int*> kromosomAs, kromosomBs, kromosomOffs;
    vector<int> odAs, odBs, mutAs, mutBs;
    vector<bool> isMuts;
    int chosenParentSize=parentsIdx.size();
    int offSize = chosenParentSize*2;
    for(int p1=0;p1<chosenParentSize;p1++){
      int pIdx1 = parentsIdx[p1];
      int* kromosomP1 = population[pIdx1]->kromosom;
      for(int p2=p1+1;p2<chosenParentSize;p2++){
        int pIdx2 = parentsIdx[p2];
        int* kromosomP2 = population[pIdx2]->kromosom;
        for(int ofIdx=0;ofIdx<2;ofIdx++){
          if (ofIdx==0){
            kromosomAs.push_back(kromosomP1);
            kromosomBs.push_back(kromosomP2);
          } else{
            kromosomAs.push_back(kromosomP2);
            kromosomBs.push_back(kromosomP1);
          }

          int* kromosomOff = create1DArrayInt(config->nCust);
          kromosomOffs.push_back(kromosomOff);

          int a=rand()%config->nCust;
          int b=rand()%config->nCust;
          if (a>b){
            int c=a;
            a=b;
            b=c;
          }
          odAs.push_back(a);
          odBs.push_back(b);

          a=rand()%config->nCust;
          b=rand()%config->nCust;
          if (a>b){
            int c=a;
            a=b;
            b=c;
          }
          mutAs.push_back(a);
          mutBs.push_back(b);

          double r=rand()/RAND_MAX;
          isMuts.push_back(r<config->pm);

        }
      }
    }

    /*
      prepare data in cuda device memory
    */
    vector<int*> d_kromosomAs, d_kromosomBs, d_kromosomOffs;
    vector<int*> d_odAs, d_odBs, d_mutAs, d_mutBs;
    vector<cudaStream_t> cudaStreams;

    for(int i=0;i<offSize;i++){
      /*
        with cuda streams, malloc and memcpy is asynchronous, hopefully faster
      */
      cudaStream_t stream;
      cudaStreamCreate(&stream);
      cudaStreams.push_back(stream);


      int *d_kromosomP1, *d_kromosomP2, *d_kromosomOff;
      d_kromosomAs.push_back(d_kromosomP1);
      d_kromosomBs.push_back(d_kromosomP2);
      d_kromosomOffs.push_back(d_kromosomOff);

      int *d_odA, *d_odB, *d_mutA, *d_mutB;
      d_odAs.push_back(d_odA);
      d_odBs.push_back(d_odB);
      d_mutAs.push_back(d_mutA);
      d_mutBs.push_back(d_mutB);

      cudaMalloc(&d_kromosomAs[i], config->nCust*sizeof(int));
      cudaMalloc(&d_kromosomBs[i], config->nCust*sizeof(int));
      cudaMalloc(&d_kromosomOffs[i], config->nCust*sizeof(int));
      cudaMalloc(&d_odAs[i], sizeof(int));
      cudaMalloc(&d_odBs[i], sizeof(int));
      cudaMalloc(&d_mutAs[i], sizeof(int));
      cudaMalloc(&d_mutBs[i], sizeof(int));

      cudaMemcpyAsync(d_kromosomAs[i], kromosomAs[i], config->nCust*sizeof(int), cudaMemcpyHostToDevice);
      cudaMemcpyAsync(d_kromosomBs[i], kromosomBs[i], config->nCust*sizeof(int), cudaMemcpyHostToDevice);
      cudaMemcpyAsync(d_odAs[i], &odAs[i], sizeof(int), cudaMemcpyHostToDevice, cudaStreams[i]);
      cudaMemcpyAsync(d_odBs[i], &odBs[i], sizeof(int), cudaMemcpyHostToDevice, cudaStreams[i]);
      cudaMemcpyAsync(d_mutAs[i], &mutAs[i], sizeof(int), cudaMemcpyHostToDevice, cudaStreams[i]);
      cudaMemcpyAsync(d_mutBs[i], &mutBs[i], sizeof(int), cudaMemcpyHostToDevice, cudaStreams[i]);
    }

    for(int i=0;i<offSize;i++){
      cudaStreamSynchronize(cudaStreams[i]);
      orderCrossover<<< 1,1,0,cudaStreams[i] >>>(d_nCust, d_kromosomAs[i], d_kromosomBs[i], d_kromosomOffs[i],
        d_odAs[i], d_odBs[i]);
    }

    for(int i=0;i<offSize;i++){
      if(isMuts[i]){
        cudaStreamSynchronize(cudaStreams[i]);
        int threadCount = mutBs[i]-mutAs[i]+1;
        rsMutationPar<<< 1,threadCount,0,cudaStreams[i] >>>(d_kromosomOffs[i],
          d_mutAs[i], d_mutBs[i]);
      }
    }

    cudaDeviceSynchronize();

    for(int i=0;i<offSize;i++){
      cudaMemcpy(kromosomOffs[i], d_kromosomOffs[i], config->nCust*sizeof(int), cudaMemcpyDeviceToHost);
      Individu *newIdv = new Individu;
      newIdv->kromosom = kromosomOffs[i];
      decodeKromosom(config, newIdv->kromosom, orderData, &newIdv->routeSet);
      calculateFitness(newIdv);
      population.push_back(newIdv);
    }

    /*
      freeing cuda memory
    */
    for(int i=0; i<offSize; i++){
      cudaFree(d_kromosomAs[i]);
      cudaFree(d_kromosomBs[i]);
      cudaFree(d_kromosomOffs[i]);
      cudaFree(d_odAs[i]);
      cudaFree(d_odBs[i]);
      cudaFree(d_mutAs[i]);
      cudaFree(d_mutBs[i]);
      cudaStreamDestroy(cudaStreams[i]);
    }
    vector<int*>().swap(d_kromosomAs);
    vector<int*>().swap(d_kromosomBs);
    vector<int*>().swap(d_kromosomOffs);
    vector<int*>().swap(d_odAs);
    vector<int*>().swap(d_odBs);
    vector<int*>().swap(d_mutAs);
    vector<int*>().swap(d_mutBs);


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
  vector<Individu*>().swap(population);

}
