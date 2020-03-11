#include<algorithm>
#include<iostream>

#include"ga_lib.h"
#include"helper_lib.h"
using namespace std;

int main(int argc, char **argv){
  ios_base::sync_with_stdio(false);

  char *configFileName = argv[1];
  Config *config = readConfig(configFileName);
  OrderData *orderData = readOrderData(config);

  /*
    Initializing Population of N individu
    initialize by random shuffes and greedy
    but the greedy is not implemented yet :D
    after initialization, evaluate and then sort by fitnes value
  */
  vector<Individu*> population;
  for(int i=0;i<config->N;i++){
    Individu* newIdv;
    if (i%2==0){
      newIdv = initIndividuRandom(config->nCust);
    } else {
      newIdv = initIndividuGreedy(config, orderData);
    }
    population.push_back(newIdv);
    population[i]->routeSet = decodeKromosom(config, population[i]->kromosom, orderData);
    calculateFitness(population[i]);
  }
  sort(population.begin(), population.end(), cmpIndividuFitness);

  /*
    Start the GA
    for MaxIter
    or until the difference between
    the difference between current bestFitness and last bestFitness
    is less than the threshold
    for 100 consecutive iterations
  */
  int sameFitnessCount = 0;
  Individu bestIndividu = *population[0];
  double bestFitness = bestIndividu.fitnessValue;
  cout<<0<<" "<<bestFitness<<" "<<bestIndividu.routeCount<<" "<<bestIndividu.totalDist<<"\n";
  for (int iter=0;iter<config->maxIter && sameFitnessCount<100;iter++){

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
      Individu* parent1 = population[pIdx1];
      for(int p2=p1+1;p2<chosenParentSize;p2++){
        int pIdx2 = parentsIdx[p2];
        Individu* parent2 = population[pIdx2];
        pair<Individu*,Individu*> parentPair = make_pair(parent1,parent2);

        pair<Individu*,Individu*> offs = orderCrossover(config, parentPair);

        /*
          mutating offspring
        */
        double r=rand()/RAND_MAX;
        if (r<config->pm){
          rsMutation(config, offs.first);
        }

        r=rand()/RAND_MAX;
        if (r<config->pm){
          rsMutation(config, offs.second);
        }
        offs.first->routeSet = decodeKromosom(config, offs.first->kromosom, orderData);
        offs.second->routeSet = decodeKromosom(config, offs.second->kromosom, orderData);
        calculateFitness(offs.first);
        calculateFitness(offs.second);
        population.push_back(offs.first);
        population.push_back(offs.second);
      }
    }

    /*
      selection NSGA2
    */
    vector<Individu*> newPopulation = selectionNSGA2(config,&population);
    population = newPopulation;


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
    cout<<iter+1<<" "<<bestFitness<<" "<<bestIndividu.routeCount<<" "<<bestIndividu.totalDist<<"\n";
  }
}
