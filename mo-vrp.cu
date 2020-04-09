#include<algorithm>
#include<iomanip>
#include<iostream>
#include<omp.h>
#include<time.h>
#include<thrust/host_vector.h>
#include<thrust/execution_policy.h>

#include"ga_lib.h"
#include"helper_lib.h"
using namespace std;

struct GenRand{
  int nCust;

  GenRand(int _nCust): nCust(_nCust){}

  int operator() (){
    return rand()%nCust;
  }
};

struct GenIsMut{
  double pm;

  GenIsMut(double _pm): pm(_pm) {}

  bool operator() (){
    double r = (double)rand()/(double)RAND_MAX;
    return r<=pm;
  }
};

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
    preparing memory
  */
  int maxOffSize = config->NP * (config->NP+1);
  thrust::host_vector<int> odAs(maxOffSize), odBs(maxOffSize);
  thrust::host_vector<int> mutAs(maxOffSize), mutBs(maxOffSize);
  thrust::host_vector<bool> isMuts(maxOffSize);
  thrust::host_vector<int*> kromosomAs(maxOffSize), kromosomBs(maxOffSize);
  thrust::host_vector<int*> kromosomOffs(maxOffSize);

  int *p_odAs = thrust::raw_pointer_cast(odAs.data());
  int *p_odBs = thrust::raw_pointer_cast(odBs.data());
  int *p_mutAs = thrust::raw_pointer_cast(mutAs.data());
  int *p_mutBs = thrust::raw_pointer_cast(mutBs.data());
  bool *p_isMuts = thrust::raw_pointer_cast(isMuts.data());
  int **p_kromosomAs = thrust::raw_pointer_cast(kromosomAs.data());
  int **p_kromosomBs = thrust::raw_pointer_cast(kromosomBs.data());
  int **p_kromosomOffs = thrust::raw_pointer_cast(kromosomOffs.data());

  // for(int i=0;i<maxOffSize;i++){
    // kromosomAs[i] = (int*) malloc(config->nCust*sizeof(int));
    // kromosomBs[i] = (int*) malloc(config->nCust*sizeof(int));
    // kromosomOffs[i] = (int*) malloc(config->nCust*sizeof(int));
  // }

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
    int chosenParentSize=parentsIdx.size();
    int offSize = chosenParentSize * (chosenParentSize - 1);

    /*
      randomizing crossover points
    */
    thrust::generate(thrust::host, odAs.begin(), odAs.begin()+offSize, GenRand(config->nCust));
    thrust::generate(thrust::host, odBs.begin(), odBs.begin()+offSize, GenRand(config->nCust));
    thrust::generate(thrust::host, mutAs.begin(), mutAs.begin()+offSize, GenRand(config->nCust));
    thrust::generate(thrust::host, mutBs.begin(), mutBs.begin()+offSize, GenRand(config->nCust));
    thrust::generate(thrust::host, isMuts.begin(), isMuts.begin()+offSize, GenIsMut(config->pm));

    int ofIdx=0;
    for(int p1=0;p1<chosenParentSize;p1++){
      int pIdx1 = parentsIdx[p1];
      int* kromosomP1 = population[pIdx1]->kromosom;
      for(int p2=p1+1;p2<chosenParentSize;p2++, ofIdx+=2){
        int pIdx2 = parentsIdx[p2];
        int* kromosomP2 = population[pIdx2]->kromosom;

        int* kromosomOff1 = create1DArrayInt(config->nCust);
        int* kromosomOff2 = create1DArrayInt(config->nCust);
        kromosomAs[ofIdx]=kromosomP1;
        kromosomBs[ofIdx]=kromosomP2;
        kromosomAs[ofIdx+1]=kromosomP2;
        kromosomBs[ofIdx+1]=kromosomP1;
        kromosomOffs[ofIdx]=kromosomOff1;
        kromosomOffs[ofIdx+1]=kromosomOff2;
      }
    }

    crossoverMutation(config->nCust, offSize, p_kromosomAs, p_kromosomBs,
      p_kromosomOffs, p_odAs, p_odBs, p_isMuts, p_mutAs, p_mutBs);

    for(int i=0;i<offSize;i++){
      Individu* offspring = new Individu();
      offspring->kromosom = kromosomOffs[i];
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
