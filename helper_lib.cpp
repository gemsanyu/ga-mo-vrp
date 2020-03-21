#include<cmath>
#include<fstream>
#include<iostream>
#include<math.h>
#include<sstream>
#include<stdlib.h>
#include<string>
#include"helper_lib.h"

using namespace std;

bool* create1DArrayBool(int sizeX){
  // bool* array = (bool*) malloc(sizeX * sizeof(bool));
  bool* array = new bool[sizeX];
  for (int i=0;i<sizeX;i++){
    array[i]=false;
  }
  return array;
}

int* create1DArrayInt(int sizeX){
    // int*  array =(int*) malloc(sizeX * sizeof(int));
    int* array = new int[sizeX];
    return array;
}

double euclideanDistance(Coordinate a, Coordinate b){
  return std::sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

Config* readConfig(string configFileName){
  ifstream configFile(configFileName);

  Config* config = new Config;
  string dataName;
  configFile >> config->maxIter >> config->threshold >>  config->N;
  configFile >> config->NP >> config->pc >> config->pm;
  configFile.close();
  return config;
}
