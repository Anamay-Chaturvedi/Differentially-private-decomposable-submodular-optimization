#include "../initialise.hpp"

class randomSet : public initialise{
public:

  randomSet  (int givenNumOfLocs1, int givenNumOfLocs2, int givenNumOfData,
                  float givenEps, float givenDelta, int givenRank1,
                  int givenRank2, vector<float> box,
                  vector<float> givenLocLat,
                  vector<float> givenLocLon);

  float algoWrapper();

};

randomSet::randomSet(int givenNumOfLocs1,
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
}

float randomSet :: algoWrapper(){

  float score[numOfLocs1 + numOfLocs2];
  float chosenSetValue = 0;

  set <int, less <int> > locations, chosenSet, nextSet;
  set <int, less <int> >::iterator itr;
  vector <float> probs(numOfLocs1 + numOfLocs2);
  int chosenIndex;

  int check1 = 0, check2 = 0;

  for(int i=0;i<numOfLocs1 + numOfLocs2;i++)
    locations.insert(i);

  if (!(rank2)){
  for(int i=0;i<rank1+rank2;i++){
    for(int j=0;j<numOfLocs1 + numOfLocs2;j++)
      probs[j] = 0;
    for(itr=locations.begin(); itr!=locations.end(); itr++)
      probs[*itr] = 1;
    std::discrete_distribution<> distributionScore(probs.begin(), probs.end());
    chosenIndex = distributionScore(generator);
    chosenSet.insert(chosenIndex);
    locations.erase(chosenIndex);

  }
  }
  else{
    for(int i=0;i<rank1+rank2;i++){
      for(int j=0;j<numOfLocs1 + numOfLocs2;j++)
        probs[j] = 0;
      for(itr=locations.begin(); itr!=locations.end(); itr++)
        probs[*itr] = 1;
      std::discrete_distribution<> distributionScore(probs.begin(), probs.end());
      chosenIndex = distributionScore(generator);
      chosenSet.insert(chosenIndex);
      locations.erase(chosenIndex);

      if (chosenIndex<numOfLocs1)
        check1++;
      else
        check2++;

      chosenSet.insert(chosenIndex);
      locations.erase(chosenIndex);

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
