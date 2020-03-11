#include<algorithm>
#include<iostream>
#include<limits>

#include"ga_lib.h"

/*
  generate pareto front layers
  based on dominationcount
*/
vector<vector<Individu*>> getParetoFronts(vector<Individu*>* population){
  int populationSize=population->size();
  for(int i=0;i<populationSize;i++){
    population->at(i)->dominatedCount=0;
    population->at(i)->dominatedIndividuVec.clear();
  }

  /*
    count domination
    and take not which dominated which
  */
  for(int i=0;i<populationSize;i++){
    for(int j=0;j<populationSize;j++){
      if (i==j){
        continue;
      }
      if (isDominate(population->at(i),population->at(j))){
        population->at(j)->dominatedCount++;
        population->at(i)->dominatedIndividuVec.push_back(population->at(j));
      }
    }
  }
  
  /*
    start from the pareto optimal front
    which consists of solutions with dominatedCount = 0
  */
  vector<Individu*>* paretoFront = new vector<Individu*>;
  for (int i=populationSize-1;i>=0;i--){
    if(population->at(i)->dominatedCount==0){
      paretoFront->push_back(population->at(i));
      population->erase(population->begin()+i);
    }
  }

  vector<vector<Individu*>> paretoFronts;

  /*
    remove the pareto optimal front from the population
    and get the next pareto optimal front from the population
    repeat
  */
  while(!population->empty()){
    paretoFronts.push_back(*paretoFront);
    vector<Individu*>* lastParetoFront = paretoFront;
    paretoFront = new vector<Individu*>;
    for(int i=0;i<lastParetoFront->size();i++){
      for(int j=0;j<lastParetoFront->at(i)->dominatedIndividuVec.size();j++){
        lastParetoFront->at(i)->dominatedIndividuVec[j]->dominatedCount--;
      }
    }

    populationSize=population->size();
    for (int i=populationSize-1;i>=0;i--){
      if(population->at(i)->dominatedCount==0){
        paretoFront->push_back(population->at(i));
        population->erase(population->begin()+i);
      }
    }
  }
  return paretoFronts;
}

/*
  select N individus from current population
  which have old population + offsprings
  by NSGAII
*/
vector<Individu*> selectionNSGA2(Config *config, vector<Individu*>* population){
  vector<Individu*> newPopulation;
  vector<vector<Individu*>> paretoFronts = getParetoFronts(population);

  bool newPopulationFull=false;
  for(int i=0;i<paretoFronts.size() && !newPopulationFull;i++){
    sortCrowdingDistance(paretoFronts[i]);
    for(int j=0;j<paretoFronts[i].size();j++){
      if(newPopulation.size()==config->N){
        newPopulationFull;
        break;
      }
      newPopulation.push_back(paretoFronts[i][j]);
    }
  }

  return newPopulation;
}

/*
  sort sub population p'
  which |p'| + |new_pop| > config.N
  so we can get the top by crowding distance
  and only add those top individus to the new_pop
*/
void sortCrowdingDistance(vector<Individu*> population){
  int populationSize = population.size();
  for(int i=0;i<populationSize;i++){
    population[i]->crowdingDistance=0;
  }

  /*
    find span (max-min) of each objective
  */
  double spanTotalDist=0.00000001, spanRouteCount=0.00000001;
  double maxTotDist=numeric_limits<double>::min(), minTotDist=numeric_limits<double>::max();
  double maxRouteCount=numeric_limits<double>::min(), minRouteCount=numeric_limits<double>::max();
  for(int i=0;i<populationSize;i++){
    if (population[i]->totalDist > maxTotDist){
      maxTotDist = population[i]->totalDist;
    }
    if (population[i]->totalDist < minTotDist){
      minTotDist = population[i]->totalDist;
    }
    if (population[i]->routeCount > maxRouteCount){
      maxRouteCount = population[i]->routeCount;
    }
    if (population[i]->routeCount < minRouteCount){
      minRouteCount = population[i]->routeCount;
    }
  }
  spanTotalDist += (maxTotDist-minTotDist);
  spanRouteCount += (maxRouteCount-minRouteCount);

  /*
    sort by crowdingDistance per objective function
  */
  sort(population.begin(), population.end(), cmpIndividuTotalDist);
  population[0]->crowdingDistance=numeric_limits<double>::max();
  population[populationSize-1]->crowdingDistance=numeric_limits<double>::max();
  for (int i=1;i<populationSize-1;i++){
    population[i]->crowdingDistance += (population[i+1]->totalDist-population[i-1]->totalDist)/spanTotalDist;
  }

  sort(population.begin(), population.end(), cmpIndividuRouteCount);
  population[0]->crowdingDistance=numeric_limits<double>::max();
  population[populationSize-1]->crowdingDistance=numeric_limits<double>::max();
  for (int i=1;i<populationSize-1;i++){
    population[i]->crowdingDistance += (population[i+1]->routeCount-population[i-1]->routeCount)/spanRouteCount;
  }

  sort(population.begin(), population.end(), cmpIndividuCrowdingDistance);
}
