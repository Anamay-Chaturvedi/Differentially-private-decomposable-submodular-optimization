#ifndef DPO_GREEDY_H
#define DPO_GREEDY_H

#include<iostream>
#include<fstream>
#include<functional>
#include<random>
#include<set>
#include<cmath>
#include<chrono>
#include "../initialise.hpp"

using namespace std;

class greedy : public initialise{
public:

  greedy  (int givenNumOfLocs1,
    int givenNumOfLocs2,
    int givenNumOfData,
    float givenEps,
    float givenDelta,
    int givenRank1,
    int givenRank2,
    vector<float> box,
    vector<float> givenLocLat,
    vector<float> givenLocLon);

  float algoWrapper();

  void timeChecker();

};

greedy::greedy  (int givenNumOfLocs1,
  int givenNumOfLocs2,
  int givenNumOfData,
  float givenEps,
  float givenDelta,
  int givenRank1,
  int givenRank2,
  vector<float> box,
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
  cout<<"Greedy initialised successfully"<<endl;
}

float greedy::algoWrapper(){

  set <int, less <int> > locations, current, prospective;
  set <int, less <int> >::iterator itr;

  float max;
  int argMax;
  float marginalGain;
  int check1 = 0, check2 = 0;

  for(int i = 0; i<numOfLocs1 + numOfLocs2; i++)
    locations.insert(i);

  //Pick the best element (rank1+rank2) many times.
  for(int i = 0; i<(rank1+rank2); i++){
    max = 0;
    //Compute the marginal benefits and keep track of the max
    for(auto itr:locations){
      prospective = current;
      prospective.insert(itr);
      marginalGain = funValue(prospective) - funValue(current);
      if (max<marginalGain){
        max = marginalGain;
        argMax = itr;
      }
    }
    if (argMax<numOfLocs1)
      check1++;
    else
      check2++;

    current.insert(argMax);
    locations.erase(argMax);

    if (check1==rank1){
      for(int j = 0; j<numOfLocs1; j++)
        locations.erase(j);
    }
    if (check2==rank2){
      for(int j = numOfLocs1; j<numOfLocs1 + numOfLocs2; j++)
        locations.erase(j);
      }

  }

  return(funValue(current));

}

#endif  // DPO_CG_H
