#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <random>
#include <vector>

#include "helper_lib.h"

void decodeKromosom(thrust::device_vector<int> kromosom, Data data,
  Config config, int &routeCount, double &totalDist){

    totalDist=0;
    routeCount=1;

    double maxDist = config.maxDist;
    int maxCap = config.maxCap;
    int lastCustIdx=-1;

    for(int i=0;i<config.nCust;i++){
      int nextCustIdx = kromosom[i];
      int orderSize = data.customers.orderSize[nextCustIdx];
      double nextCustX = data.customers.x[nextCustIdx];
      double nextCustY = data.customers.y[nextCustIdx];
      double dist;
      if (lastCustIdx==-1){
        dist = data.distancesToDepot[nextCustIdx];
      } else {
        dist = data.distancesToCust[lastCustIdx*config.nCust+nextCustIdx];
      }
      double distToDepot = data.distancesToDepot[nextCustIdx];

      /*
        check if violation of max dist || max cap happens
        then go back to depot first
      */
      if (maxDist<dist+distToDepot || maxCap<orderSize){
        totalDist += data.distancesToDepot[lastCustIdx];
        routeCount++;
        maxCap = config.maxCap;
        maxDist = config.maxDist;

        dist = data.distancesToDepot[nextCustIdx];
      }

      maxDist-=dist;
      maxCap-=orderSize;
      totalDist+=dist;
      lastCustIdx=nextCustIdx;
    }
}
