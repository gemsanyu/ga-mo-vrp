#include<algorithm>
#include<fstream>
#include<iostream>
#include<limits>
#include<omp.h>
#include<random>
#include<sstream>
#include<stdlib.h>
#include<string>
#include<utility>

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

bool cmpIndividuFitness(Individu* a, Individu* b){
  return a->fitnessValue > b->fitnessValue;
}

bool cmpIndividuCrowdingDistance(Individu* a, Individu* b){
  return a->crowdingDistance > b->crowdingDistance;
}

bool cmpIndividuTotalDist(Individu* a, Individu* b){
  return a->totalDist < b->totalDist;
}

bool cmpIndividuRouteCount(Individu* a, Individu* b){
  return a->routeCount < b->routeCount;
}

Individu* create1DArrayIndividu(int size){
  Individu *array = new Individu[size];
  return array;
}

Customer* create1DArrayCustomer(int size){
  Customer *array = new Customer[size];
  return array;
}

void decodeKromosom(Config *config, int *kromosom, OrderData *orderData, RouteSet* routeSet){
  std::vector<int> route;
  double totalDist=0;
  int totalOrder=0;
  Coordinate lastCoord=orderData->depot;
  for (int k=0;k<config->nCust;k++){
    int custID = kromosom[k];
    double dist = euclideanDistance(lastCoord, orderData->customerData[custID].coordinate);
    double distToDepot = euclideanDistance(orderData->customerData[custID].coordinate, orderData->depot);
    int odSize = orderData->customerData[custID].orderSize;
    if ((totalOrder+odSize>config->maxCap) || (totalDist+dist+distToDepot>config->maxDist)){
      routeSet->routes.push_back(route);
      totalDist += euclideanDistance(lastCoord, orderData->depot);
      routeSet->distances.push_back(totalDist);
      route.clear();
      totalDist=0;
      totalOrder=0;
      lastCoord=orderData->depot;
      dist = euclideanDistance(lastCoord, orderData->customerData[custID].coordinate);
    }
    route.push_back(custID);
    totalDist += dist;
    totalOrder += orderData->customerData[custID].orderSize;
    lastCoord = orderData->customerData[custID].coordinate;
  }
  routeSet->routes.push_back(route);
  totalDist += euclideanDistance(lastCoord, orderData->depot);
  routeSet->distances.push_back(totalDist);
}

int* encodeRouteSet(Config *config, RouteSet *routeSet){
  int* kromosom = create1DArrayInt(config->nCust);
  int custCount=0;
  int routeCount = routeSet->routes.size();
  for(int rIdx=0;rIdx<routeCount;rIdx++){
    vector<int> route = routeSet->routes[rIdx];
    int routeLength=route.size();
    for(int cIdx=0;cIdx<routeLength;cIdx++){
      kromosom[custCount]=route[cIdx];
      custCount++;
    }
  }
  return kromosom;
}

void initIndividuRandom(int nCust, int* kromosom){
  for(int i=0;i<nCust;i++){
    kromosom[i]=i;
  }
  random_shuffle(kromosom, kromosom+nCust);
}

/*
  find nearest feasible customer
  feasible if cap still can contain ordersize
  and if the distance to the customer and back to depot + current dist <= config.maxdist
*/
int findNearestCustIdx(Config* config, OrderData* orderData, vector<int>* custsIdx, Coordinate* currentCoord, int currentCap, double currentDist){
  double maxDist = config->maxDist - currentDist;
  int maxCap = config->maxCap - currentCap;
  Coordinate depotCoord = orderData->depot;

  double closestDist = numeric_limits<double>::max();
  int closestIdx = -1;
  for(int i=0;i<custsIdx->size();i++){
    int custIdx = custsIdx->at(i);

    /*
      checking capacity constraint
    */
    int odSize = orderData->customerData[custIdx].orderSize;
    if(odSize>maxCap){
      continue;
    }

    /*
      checking distance and dist to depot constraint
    */
    Coordinate custCoord = orderData->customerData[custIdx].coordinate;
    double dist = euclideanDistance(*currentCoord, custCoord);
    double distToDepot = euclideanDistance(custCoord, depotCoord);
    if (dist+distToDepot>maxDist){
      continue;
    }

    if (dist<closestDist){
      closestIdx = i;
      closestDist = dist;
    }
  }

  return closestIdx;
}

void initIndividuGreedy(Config* config, OrderData* orderData, int* kromosom){
  vector<int> custsIdx;
  for(int i=0;i<config->nCust;i++){
    custsIdx.push_back(i);
  }

  int servedCustCount = 0;
  double totalDist=0;
  int totalOrder=0;
  Coordinate lastCoord=orderData->depot;
  while(!custsIdx.empty()){
    //randomize the first customer
    int closestIdx;
    if (servedCustCount == 0){
      closestIdx = rand()%config->nCust;
    } else {
      closestIdx = findNearestCustIdx(config, orderData, &custsIdx, &lastCoord,
        totalOrder, totalDist);
    }

    /*
      if no feasible customer left to serve
      reset everything
    */
    if (closestIdx == -1){
      totalDist = 0;
      totalOrder = 0;
      lastCoord=orderData->depot;
      continue;
    }

    int closestCustIdx=custsIdx[closestIdx];
    totalOrder+= orderData->customerData[closestCustIdx].orderSize;
    Coordinate custCoord = orderData->customerData[closestCustIdx].coordinate;
    double dist = euclideanDistance(lastCoord, custCoord);
    totalDist += dist;

    lastCoord = custCoord;
    kromosom[servedCustCount]=closestCustIdx;
    servedCustCount++;
    custsIdx.erase(custsIdx.begin()+closestIdx);
  }
}

bool isDominate(Individu* idvA, Individu* idvB){
  return (idvA->totalDist <= idvB->totalDist) &&
  (idvA->routeCount<=idvB->routeCount) &&
  ((idvA->totalDist<idvB->totalDist) || (idvA->routeCount<idvB->routeCount));
}

void readOrderData(Config *config, OrderData *odData){
  ifstream dataFile(config->fileName);
  dataFile >> config->nCust >> config->maxDist >> config->maxCap;
  dataFile >> odData->depot.x >> odData->depot.y;
  odData->customerData = create1DArrayCustomer(config->nCust);
  for (int c=0;c<config->nCust;c++){
    dataFile >> odData->customerData[c].coordinate.x;
    dataFile >> odData->customerData[c].coordinate.y;
    dataFile >> odData->customerData[c].orderSize;
  }
  dataFile.close();
}
