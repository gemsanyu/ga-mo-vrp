#ifndef GA_LIB_H
#define GA_LIB_H
#include<string>
#include<utility>
#include<vector>

#include"helper_lib.h"
using namespace std;

struct RouteSet{
  vector<vector<int> > routes;
  vector<double> distances;
};

struct Individu{
  int *kromosom;
  RouteSet routeSet;
  double totalDist;
  int routeCount;
  double fitnessValue;
  double crowdingDistance;
  int dominatedCount;
  vector<Individu*> dominatedIndividuVec;
};

struct Customer{
  Coordinate coordinate;
  int orderSize;
};

struct OrderData{
  Customer *customerData;
  Coordinate depot;
};

void calculateFitness(Individu *individu);
bool cmpIndividuFitness(Individu* a, Individu* b);
bool cmpIndividuCrowdingDistance(Individu* a, Individu* b);
bool cmpIndividuTotalDist(Individu* a, Individu* b);
bool cmpIndividuRouteCount(Individu* a, Individu* b);
Individu* create1DArrayIndividu(int size);
Customer* create1DArrayCustomer(int size);
void decodeKromosom(Config *config, int *kromosom, OrderData *orderData, RouteSet* routeSet);
int* encodeRouteSet(Config *config, RouteSet *routeSet);
void initIndividuRandom(int nCust, int* kromosom);
void initIndividuGreedy(Config* config, OrderData* orderData, int* kromosom);
bool isDominate(Individu* idvA, Individu* idvB);
void readOrderData(Config *config, OrderData *orderData);
vector<Individu*> selectionNSGA2(Config *config, vector<Individu*>* population);
void sortCrowdingDistance(vector<Individu*>population);

#endif
