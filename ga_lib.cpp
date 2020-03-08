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

pair<Individu,Individu> orderCrossover(pair<Individu,Individu> parents){

}

void rsmMutation(Individu *Individu){
  
}

Individu orderCrossover_(Individu parentA, Individu parentB){

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
