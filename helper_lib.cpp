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
  bool* array = (bool*) malloc(sizeX * sizeof(bool));
  for (int i=0;i<sizeX;i++){
    array[i]=false;
  }
  return array;
}

int* create1DArrayInt(int sizeX){
    int*  array =(int*) malloc(sizeX * sizeof(int));
    return array;
}

double euclideanDistance(Coordinate a, Coordinate b){
  return std::sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

Config readConfig(string configFileName){
  ifstream configFile(configFileName);

  Config config;
  string dataName;
  configFile >> dataName;
  config.fileName = "data/"+dataName;
  configFile >> config.nCust >> config.maxDist >> config.maxCap >> config.maxIter;
  configFile >> config.threshold >>  config.N >> config.NP >> config.pc >> config.pm;
  cout << "Using " << dataName <<"\n";
  cout << "Num of Cust : " << config.nCust <<"\n";
  cout << "Max Dist : " << config.maxDist <<"\n";
  cout << "Max Cap : " << config.maxCap <<"\n";
  cout << "Max Iteration : " << config.maxIter  <<"\n";
  cout << "convergence threshold : " << config.threshold  <<"\n";
  cout << "N : " << config.N  <<"\n";
  cout << "NP : " << config.NP<<"\n";
  cout << "pc : " << config.pc<<"\n";
  cout << "pm : " << config.pm<<"\n";
  configFile.close();
  return config;
}

/*
CONFIG SMALL-DISTRIBUTED-1
Using small-distributed
Num of Cust : 50
Max Dist : 200
Max Cap : 30
Max Iteration : 1000000
convergence threshold : 1e-09
N : 1
NP : 10
pc : 0.8
pm : 0.03
*/
