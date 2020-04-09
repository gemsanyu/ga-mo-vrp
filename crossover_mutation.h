#ifndef CM_H
#define CM_H
#include<string>
#include<utility>
#include<vector>
#include<thrust/host_vector.h>

#include"helper_lib.h"
#include"ga_lib.h"
using namespace std;

void spinRouletteWheel(vector<Individu*>* population, int spinCount,
  thrust::host_vector<double> &probs, thrust::host_vector<int> &parentIdxs);
void crossoverMutation(int nCust, int offSize, int **kromosomAs, int **kromosomBs,
  int **kromosomOffs, int *odAs, int *odBs, bool *isMuts, int *mutAs, int *mutBs);

#endif
