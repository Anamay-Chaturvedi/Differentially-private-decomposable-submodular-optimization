#ifndef DPO_DIFFPRIVATEGREEDY_H
#define DPO_DIFFPRIVATEGREEDY_H

#include<iostream>
#include<fstream>
#include<functional>
#include<random>
#include<set>
#include<cmath>
#include<chrono>
#include "../initialise.hpp"

using namespace std;


class diffPrivateGreedy : public initialise{
public:
  float eps_0;

  diffPrivateGreedy  (int givenNumOfLocs1, int givenNumOfLocs2, int givenNumOfData,
                  float givenEps, float givenDelta, int givenRank1,
                  int givenRank2, vector<float> box, bool guptaPrivacy,
                  vector<float> givenLocLat,
                  vector<float> givenLocLon);

  float algoWrapper();


};

diffPrivateGreedy :: diffPrivateGreedy  (int givenNumOfLocs1,
  int givenNumOfLocs2,
  int givenNumOfData,
  float givenEps,
  float givenDelta,
  int givenRank1,
  int givenRank2,
  vector<float> box,
  bool guptaPrivacy,
  vector<float> givenLocLat,
  vector<float> givenLocLon)
  : initialise (givenNumOfLocs1,
    givenNumOfLocs2,
    givenNumOfData,
    givenEps,
    givenDelta,
    givenRank1,
    givenRank2,
    box,
    givenLocLat,
    givenLocLon) {

  eps=givenEps;
  delta=givenDelta;
  int totalRank = givenRank1 + givenRank2;
  cout<<"Total rank is "<<totalRank;
  if (!(guptaPrivacy)){
    float eps_00 = eps/totalRank;
    float eps_01 = sqrt((2*log(1.0/delta)/totalRank) + (eps/totalRank)) - sqrt((2*log(1.0/delta)/totalRank));
    if (eps_00>eps_01){
      eps_0 = eps_00;
      cout<<endl<<"Basic composition "<<endl;
    }
    else{
      eps_0 = eps_01;
      cout<<endl<<"Advanced composition "<<endl;
    }
  }
  else
    eps_0 = 2*eps/((e-1)*(3 + log(1/delta)));


  cout<<endl<<"dPG eps_0 is "<<eps_0<<endl;
}

float diffPrivateGreedy :: algoWrapper(){

  float score[numOfLocs1 + numOfLocs2];
  float chosenSetValue = 0;

  int check1 = 0;
  int check2 = 0;

  for(int i=0;i<(numOfLocs1 + numOfLocs2);i++)
    score[i] = 0.0;

  set <int, less <int> > locations, chosenSet, nextSet;
  set <int, less <int> >::iterator itr;
  vector <float> probs(numOfLocs1 + numOfLocs2);
  int chosenIndex;


  for(int i=0;i<numOfLocs1 + numOfLocs2;i++)
    locations.insert(i);

  if (!(rank2)){
    for(int i=0;i<rank1;i++){
      chosenSetValue = funValue(chosenSet);
      for(itr=locations.begin(); itr!=locations.end(); itr++){
        nextSet = chosenSet;
        nextSet.insert(*itr);
        score[*itr] = funValue(nextSet) - chosenSetValue;
        probs[*itr] = exp(eps_0*(score[*itr])/2.0);
      }
      for(itr=chosenSet.begin(); itr!=chosenSet.end(); itr++)
        probs[*itr] = 0;

      std::discrete_distribution<> distributionScore(probs.begin(), probs.end());
      chosenIndex = distributionScore(generator);
      chosenSet.insert(chosenIndex);
    }
  }

  else{
  for(int i=0;i<rank1+rank2;i++){
    chosenSetValue = funValue(chosenSet);
    for(int j=0;j<numOfLocs1 + numOfLocs2;j++)
      probs[j] = 0;
    for(itr=locations.begin(); itr!=locations.end(); itr++){
      nextSet = chosenSet;
      nextSet.insert(*itr);
      score[*itr] = funValue(nextSet) - chosenSetValue;
      probs[*itr] = exp(eps_0*(score[*itr])/2.0);
    }
    for(itr=chosenSet.begin(); itr!=chosenSet.end(); itr++)
      probs[*itr] = 0;
    std::discrete_distribution<> distributionScore(probs.begin(), probs.end());
    chosenIndex = distributionScore(generator);
    chosenSet.insert(chosenIndex);

    if (chosenIndex< numOfLocs1)
      check1++;
    else
      check2++;

    if (check1==rank1)
      for(int j = 0; j<numOfLocs1; j++)
        locations.erase(j);
    if (check2==rank2)
      for(int j = numOfLocs1; j<numOfLocs1 + numOfLocs2; j++)
        locations.erase(j);

  }

}
return(funValue(chosenSet));
}



#endif  // DPO_CG_H
