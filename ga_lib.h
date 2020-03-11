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
  RouteSet *routeSet;
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
RouteSet* decodeKromosom(Config *config, int *kromosom, OrderData *orderData);
Individu* initIndividuRandom(int nCust);
Individu* initIndividuGreedy(Config* config, OrderData* orderData);
bool isDominate(Individu* idvA, Individu* idvB);
pair<Individu*,Individu*> orderCrossover(Config* config, pair<Individu*,Individu*> parents);
void rsMutation(Config* config, Individu *individu);
OrderData* readOrderData(Config *config);
vector<Individu*> selectionNSGA2(Config *config, vector<Individu*>* population);
void sortCrowdingDistance(vector<Individu*>population);
vector<int> spinRouletteWheel(vector<Individu*>* population, int spinCount);

#endif
