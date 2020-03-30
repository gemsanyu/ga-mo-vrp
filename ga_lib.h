#ifndef GA_LIB_H
#define GA_LIB_H
#include<string>
#include<utility>
#include<vector>

#include"helper_lib.h"
using namespace std;

struct RouteSet{
  vector<vector<int>> routes;
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


#endif
