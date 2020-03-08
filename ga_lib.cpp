#include<algorithm>
#include<fstream>
#include<iostream>
#include<sstream>
#include<stdlib.h>
#include<string>
#include<utility>
#include<random>

#include"ga_lib.h"
#include"helper_lib.h"

void calculateFitness(Individu *individu){
  double totalDist = 0;
  int routeCount=individu->routeSet.routes.size();
  for (int r=0;r<routeCount;r++){
    totalDist += individu->routeSet.distances[r];
  }
  double fitnessValue = 1.0/(600.0*totalDist + 10000.0*(double)routeCount);
  individu->totalDist = totalDist;
  individu->routeCount = routeCount;
  individu->fitnessValue = fitnessValue;
}

bool cmpIndividuFitness(Individu a, Individu b){
  return a.fitnessValue < b.fitnessValue;
}

Individu* create1DArrayIndividu(int size){
  Individu *array = new Individu[size];
  return array;
}

Customer* create1DArrayCustomer(int size){
  Customer *array = new Customer[size];
  return array;
}

RouteSet decodeKromosom(Config config, int *kromosom, OrderData orderData){
  RouteSet routeSet;
  std::vector<int> route;
  double totalDist=0;
  int totalOrder=0;
  Coordinate lastCoord=orderData.depot;
  for (int k=0;k<config.nCust;k++){
    int custID = kromosom[k];
    double dist = euclideanDistance(lastCoord, orderData.customerData[custID].coordinate);
    double distToDepot = euclideanDistance(orderData.customerData[custID].coordinate, orderData.depot);
    int odSize = orderData.customerData[custID].orderSize;
    if ((totalOrder+odSize>config.maxCap) || (totalDist+dist+distToDepot>config.maxDist)){
      routeSet.routes.push_back(route);
      totalDist += euclideanDistance(lastCoord, orderData.depot);
      routeSet.distances.push_back(totalDist);
      route.clear();
      totalDist=0;
      totalOrder=0;
    }
    route.push_back(custID);
    totalDist += dist;
    totalOrder += orderData.customerData[custID].orderSize;
    lastCoord = orderData.customerData[custID].coordinate;
  }
  routeSet.routes.push_back(route);
  totalDist += euclideanDistance(lastCoord, orderData.depot);
  routeSet.distances.push_back(totalDist);
  return routeSet;
}

Individu initIndividuRandom(int nCust){
  Individu individu;
  individu.kromosom = create1DArrayInt(nCust);
  for(int i=0;i<nCust;i++){
    individu.kromosom[i]=i;
  }
  random_shuffle(individu.kromosom, individu.kromosom+nCust);
  return individu;
}

Individu orderCrossover_(Config config, Individu parentA, Individu parentB){
  /*
    First randomize segment points a and b
  */
  int a=rand()%config.nCust;
  int b=rand()%config.nCust;
  if (a>b){
    int c=a;
    a=b;
    b=c;
  }

  /*
    copy parentA segment (a,b) to offspring segment(a,b)
  */
  bool *genExistFlag = create1DArrayBool(config.nCust);
  Individu offspring;
  offspring.kromosom = create1DArrayInt(config.nCust);
  for (int c=a;c<=b;c++){
    int custID = parentA.kromosom[c];
    offspring.kromosom[c] = custID;
    genExistFlag[custID] = true;
  }

  /*
    and then add parentB's gens
    not yet contained by the offspring to offspring
  */
  int ofIdx=b+1;
  for (int genBIdx=b+1;genBIdx<a;genBIdx = (genBIdx+1)%config.nCust){
    int gen = parentB.kromosom[genBIdx];
    if (genExistFlag[gen]){
      continue;
    }
    offspring.kromosom[ofIdx]=gen;
    ofIdx = (ofIdx+1)%config.nCust;
  }

  return offspring;
}

pair<Individu,Individu> orderCrossover(Config config, pair<Individu,Individu> parents){
  pair<Individu,Individu> offs;
  offs.first = orderCrossover_(config, parents.first, parents.second);
  offs.second = orderCrossover_(config, parents.second, parents.first);
}

OrderData readOrderData(Config config){
  ifstream dataFile(config.fileName);
  OrderData odData;
  odData.customerData = create1DArrayCustomer(config.nCust);
  dataFile >> odData.depot.x >> odData.depot.y;
  for (int c=0;c<config.nCust;c++){
    dataFile >> odData.customerData[c].coordinate.x;
    dataFile >> odData.customerData[c].coordinate.y;
    dataFile >> odData.customerData[c].orderSize;
  }
  dataFile.close();
  return odData;
}
