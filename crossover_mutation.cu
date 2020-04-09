#include<thrust/execution_policy.h>
#include<thrust/host_vector.h>
#include<thrust/reduce.h>
#include<thrust/transform.h>
#include<thrust/scan.h>
using namespace std;

#include"crossover_mutation.h"
#include"ga_lib.h"
#include"helper_lib.h"

struct DivideBy{
  double total;

  DivideBy(double _total): total(_total) {}

  double operator () (double r){
    return r/total;
  }
};

// Roulette Wheel
int spinRouletteWheel_(thrust::host_vector<double> const &probs){
	double 	select; 	// Variable Random Number

	//Random Number
	select = (double) rand()/RAND_MAX;
	int c = 0;
	//code selected colom
	while(c<probs.size()){
    select-=probs[c];
		if (select <= 0){
			break;
		}
		c++;
	}
	return c;
}

void spinRouletteWheel(vector<Individu*>* population, int spinCount,
  thrust::host_vector<double> &probs, thrust::host_vector<int> &parentIdxs){

  for(int i=0;i<population->size();i++){
    probs[i]=population->at(i)->fitnessValue;
  }
  double sumProb = thrust::reduce(
    thrust::host,
    probs.begin(),
    probs.end(),
    0.
  );

  //normalize
  thrust::transform(
    thrust::host,
    probs.begin(),
    probs.end(),
    probs.begin(),
    DivideBy(sumProb)
  );

  thrust::inclusive_scan(
    thrust::host,
    probs.begin(),
    probs.end(),
    probs.begin()
  );

  for(int i=0;i<spinCount;i++){
    int res=spinRouletteWheel_(probs);
    parentIdxs[i]=res;
  }
}


void crossoverMutation(int nCust, int offSize, int **kromosomAs, int **kromosomBs,
  int **kromosomOffs, int *odAs, int *odBs, bool *isMuts, int *mutAs, int *mutBs){

  bool *genExistFlag = create1DArrayBool(nCust);
  for(int idx=0;idx<offSize;idx++){
    int odA = odAs[idx];
    int odB = odBs[idx];
    int mutA = mutAs[idx];
    int mutB = mutBs[idx];

    int *kromosomA = kromosomAs[idx];
    int *kromosomB = kromosomBs[idx];
    int *kromosomOff = kromosomOffs[idx];

    // cout<<odA<<" "<<odB<<" "<<isMut<<" "<<mutA<<" "<<mutB<<"\n";
    if (odA>odB){
      int c=odA;
      odA=odB;
      odB=c;
    }

    /*
      copy parentA segment (a,b) to offspring segment(a,b)
    */
    for(int i=0;i<nCust;i++){
      genExistFlag[i]=false;
    }

    for (int c=odA;c<=odB;c++){
      int custID = kromosomA[c];
      kromosomOff[c] = custID;
      genExistFlag[custID] = true;
    }

    /*
      and then add parentB's gens
      not yet contained by the offspring
    */
    int ofIdx=(odB+1)%nCust;
    for (int genBIdx=(odB+1)%nCust;ofIdx<odA || ofIdx>odB;genBIdx = (genBIdx+1)%nCust){
      int gen = kromosomB[genBIdx];
      if (genExistFlag[gen]){
        continue;
      }
      kromosomOff[ofIdx]=gen;

      genExistFlag[gen]=true;
      ofIdx = (ofIdx+1)%nCust;
    }

    if(!isMuts[idx]){
      continue;
    }

    /*
      swapping mutation
    */
    if(mutA>mutB){
      int c = mutA;
      mutA = mutB;
      mutB = c;
    }

    int halfCount = (mutB-mutA+1)/2;
    for(int i=0;i<halfCount;i++){
      int custID = kromosomOff[mutA+i];
      kromosomOff[mutA+i] = kromosomOff[mutB-i];
      kromosomOff[mutB-i] = custID;
    }
  }

  delete[] genExistFlag;
}
