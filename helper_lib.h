#ifndef HLP_LIB_H
#define HLP_LIB_H
#include<string>
using namespace std;

struct Config{
  string fileName;
  int nCust;
  int maxDist, maxCap;
  double threshold;
  int maxIter, N, NP;
  double pc, pm;
};


struct Coordinate{
  double x,y;
};

bool* create1DArrayBool(int sizeX);
int* create1DArrayInt(int sizeX);
double euclideanDistance(Coordinate a, Coordinate b);
Config* readConfig(string configFileName);
#endif
