#include<algorithm>
#include<iomanip>
#include<iostream>
#include<omp.h>
#include<time.h>

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

  cout << "Using Data : " << config->fileName  <<"\n";
  cout << "N Customer : " << config->nCust  <<"\n";
  cout << "Max Distance : " << config->maxDist  <<"\n";
  cout << "Max Cap : " << config->maxCap  <<"\n";
  cout << fixed << setprecision(5)<< orderData->depot.x<<" "<<orderData->depot.y<<"\n";
  cout << "Max Iteration : " << config->maxIter  <<"\n";
  cout << "convergence threshold : " << config->threshold  <<"\n";
  cout << "N : " << config->N  <<"\n";
  cout << "NP : " << config->NP<<"\n";
  cout << "pc : " << config->pc<<"\n";
  cout << "pm : " << config->pm<<"\n";

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
    int chosenParentSize=parentsIdx.size();
    for(int p1=0;p1<chosenParentSize;p1++){
      int pIdx1 = parentsIdx[p1];
      int* kromosomP1 = population[pIdx1]->kromosom;
      for(int p2=p1+1;p2<chosenParentSize;p2++){
        int pIdx2 = parentsIdx[p2];
        int* kromosomP2 = population[pIdx2]->kromosom;

        int* kromosomOff1 = create1DArrayInt(config->nCust);
        orderCrossover(config, kromosomP1, kromosomP2, kromosomOff1);

        int* kromosomOff2 = create1DArrayInt(config->nCust);
        orderCrossover(config, kromosomP2, kromosomP1, kromosomOff2);
        /*
          mutating offspring
        */
        double r=rand()/RAND_MAX;
        if (r<config->pm){
          rsMutation(config, kromosomOff1);
        }

        r=rand()/RAND_MAX;
        if (r<config->pm){
          rsMutation(config, kromosomOff2);
        }

        Individu* off1 = new Individu();
        off1->kromosom = kromosomOff1;
        decodeKromosom(config, off1->kromosom, orderData, &off1->routeSet);
        calculateFitness(off1);

        Individu* off2 = new Individu;
        off2->kromosom = kromosomOff2;
        decodeKromosom(config, off2->kromosom, orderData, &off2->routeSet);
        calculateFitness(off2);
        population.push_back(off1);
        population.push_back(off2);
      }
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

  cout<<"Obtained PF\n";
  for(int i=0;i<population.size();i++){
    cout<<fixed<<setprecision(8)<<"solution-"<<i<<": "<<population[i]->routeCount<<" "<<population[i]->totalDist<<"\n";
  }
  cout<<fixed<<setprecision(8)<<"Initial Population Generation Time: "<<initTime<<"\n";
  cout<<"GA stopped at iteration-"<<iter<<"\n";
  cout<<fixed<<setprecision(8)<<"Total Time: "<<totalTime<<"\n";
  cout<<fixed<<setprecision(8)<<"Mean Time per iteration: "<<totalTime/(double)iter<<"\n";

  /*
    freeing memory
  */
  vector<Individu*>().swap(population);
}
