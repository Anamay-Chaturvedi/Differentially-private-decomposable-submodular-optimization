#ifndef DPO_CG_H
#define DPO_CG_H

#include<iostream>
#include<fstream>
#include<functional>
#include<random>
#include<set>
#include<cmath>
#include<chrono>
#include "../initialise.hpp"

using namespace std;

class markedSamp{
public:
    int index;
    double val;

    markedSamp();
    markedSamp(int, double);
    markedSamp( const markedSamp &obj);
    ~markedSamp();
};

markedSamp::markedSamp(){}

markedSamp::markedSamp(int givenIndex, double givenVal){
  index=givenIndex;
  val=givenVal;
}

markedSamp::markedSamp( const markedSamp &obj){
  index=obj.index;
  val=obj.val;
}

markedSamp::~markedSamp(){}


class continuousGreedy : public initialise{
  public:
    float ctsScore;
    set <int, less<int> > returnSet;
    float **hashes;
    float *scores;
    float *point;
    float *prospectivePoint;
    int numOfSamples;
    int T;
    double tolerance;
    float eta;
    std::vector<float> probs;
    float probSum=0.0;
    float probMax=0.0;
    bool direct;


    std::uniform_real_distribution<double> distribution;

    continuousGreedy  (int givenNumOfLocs1, int givenNumOfLocs2, int givenNumOfData,
                    int givenNumOfSamples, int givenT, float givenEps,float givenDelta,
                    int givenRank1, int givenRank2, vector<float> box,
                    vector<float> givenLocLat,
                    vector<float> givenLocLon,
                    bool givenDirect);

    ~continuousGreedy();

    void scoreComputer();
    void scoreComputer1();
    void scoreComputer2();
    float mult(float*);

    float score(int);
    float score1(int);
    float score2(int);


    void hashRenew();

    void exponentialMechanism();
    set<int, less <int> > exponentialMechanism2();

    static bool markedSampComp(markedSamp,markedSamp);

    float algoWrapper();
    float algoWrapper2();

    std::set<int, std::less<int> > swapRounding();

    std::set<int, std::less<int> > partitionRound(float *point, int n);

    void timeChecker();

};

continuousGreedy::~continuousGreedy(){
  for(int i=0; i<numOfLocs1 + numOfLocs2; i++)
    delete[] hashes[i];
  delete[] hashes;
  delete scores;
  delete point;
  delete prospectivePoint;
}

continuousGreedy::continuousGreedy  (int givenNumOfLocs1,
  int givenNumOfLocs2,
  int givenNumOfData,
  int givenNumOfSamples,
  int givenT,
  float givenEps,
  float givenDelta,
  int givenRank1,
  int givenRank2,
  vector<float> box,
  vector<float> givenLocLat,
  vector<float> givenLocLon,
  bool givenDirect)
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
      float delta = 1.0/(2*pow(givenNumOfData,1.5));
      eps = 2*log(1 + givenEps/(2*(3 - log(delta))));

      numOfSamples = givenNumOfSamples;

      //Generate hash values to compute multilinear relaxation.
      if (!(direct)){
      float **tempHashes = new float*[numOfLocs1+numOfLocs2];

      for(int i=0;i<(numOfLocs1+numOfLocs2);i++){
        tempHashes[i] = new float[numOfSamples];
        for(int j=0;j<numOfSamples;j++)
          tempHashes[i][j] = distribution(generator);
      }
      //8
      hashes = tempHashes;
      }
      //Hash values computed.

      T = givenT;
      eta = 1.0/T;

      //10
      scores = new float[numOfLocs1 + numOfLocs2];

      //11
      point = new float[numOfLocs1 + numOfLocs2];
      prospectivePoint = new float[numOfLocs1 + numOfLocs2];

      for(int i=0;i<numOfLocs1 + numOfLocs2;i++)
        point[i] = 0;

      //12
      std::uniform_real_distribution<double> distribution1(0.0,1.0);
      distribution = distribution1;

      //Rejection sampling parameters
      tolerance = exp(eps*(sqrt((3*(rank1+rank2) + 2.5)*47 +14) ));
      cout<<"Tolerance is "<<tolerance<<endl;

      direct = givenDirect;

}

void continuousGreedy::hashRenew(){
  for(int i=0;i<(numOfLocs1 + numOfLocs2);i++)
    for(int j=0;j<numOfSamples;j++)
      hashes[i][j] = distribution(generator);
}


float continuousGreedy::score(int i){

  if (direct)
    return score2(i);
  else
    return score1(i);
}

void continuousGreedy::scoreComputer(){
  if (direct)
    scoreComputer2();
  else
    scoreComputer1();
}

float continuousGreedy::score1(int i){

  float score=0;

  for(int j=0; j<numOfSamples; j++)
    if ((hashes[i][j]>point[i])&&(hashes[i][j]<= (point[i] + eta))){
      set<int, less<int> > randomSetBig, randomSetSmall;
      for(int k=0; k<(numOfLocs1 + numOfLocs2); k++)
        if (hashes[k][j]<=point[k])
          randomSetSmall.insert(k);
      randomSetBig = randomSetSmall;
      randomSetBig.insert(i);
      score+= (funValue(randomSetBig) - funValue(randomSetSmall));
    }
  score = score/numOfSamples;
  return score;
}

void continuousGreedy::scoreComputer1(){

  //Generate a tally of which samples are needed for a coordinate's score.

  for(int i=0; i<(numOfLocs1 + numOfLocs2); i++)
    scores[i] = score1(i);

  probs.clear();
  probs.resize(numOfLocs1+numOfLocs2,0);
  probSum=0.0;
  probMax=0.0;
  for(int i=0;i<(numOfLocs1 + numOfLocs2);i++){
    probs[i] = (exp(eps*(scores[i])/2));
    probSum += probs[i];
    if(probs[i]>probMax)
      probMax=probs[i];
  }

  /*
  cout<<"score 1, number of samples "<<numOfSamples<<"\t \t";
  for(int i=0; i<(numOfLocs1 + numOfLocs2); i++){
    if (i<15)
      cout<<i<<":"<<scores[i]<<" ";
  }
  cout<<endl;
  */

}

void continuousGreedy::scoreComputer2(){

  for(int i=0;i<numOfLocs1+numOfLocs2;i++)
    scores[i] = score2(i);

  probs.clear();
  probs.resize(numOfLocs1+numOfLocs2,0);
  probSum=0.0;
  probMax=0.0;
  for(int i=0;i<(numOfLocs1 + numOfLocs2);i++){
    probs[i] = (exp(eps*(scores[i])/2));
    probSum += probs[i];
    if(probs[i]>probMax)
      probMax=probs[i];
  }
  /*
  cout<<"score 2, number of samples "<<numOfSamples<<"\t \t";
  for(int i=0; i<(numOfLocs1 + numOfLocs2); i++){
    if (i<15)
      cout<<i<<":"<<scores[i]<<" ";
  }
  cout<<endl;
  */

}

float continuousGreedy::score2(int i){
  float score=0;

  for(int j=0;j<numOfLocs1+numOfLocs2;j++)
    prospectivePoint[j] = point[j];
  prospectivePoint[i] += eta;
  score = mult(prospectivePoint) - mult(point);

  return score;
}

float continuousGreedy::mult(float* point){
  float ans=0;
  float prod=1.0;
  for(int i=0;i<numOfData;i++){
    prod=1.0;
    for(int j=0;j<numOfLocs1+numOfLocs2;j++){
      ans+=((1.0-allDistances[i][orderedLocations[i][j]])*point[orderedLocations[i][j]]*prod);
      prod = prod*(1.0 - point[orderedLocations[i][j]]);
    }
  }
  return ans;
}

void continuousGreedy::exponentialMechanism(){

    for(int r = 0;r<T; r++){
      int total = 0;
      bool chosen[numOfLocs1+numOfLocs2];
      for(int i=0;i<numOfLocs1+numOfLocs2;i++)
        chosen[i]=false;

      while(total<rank1+rank2){
        scoreComputer();
        int M = 5*(rank1+rank2-total);
        vector<markedSamp> samples (M);

        for(int i=0;i<M;i++){
          samples[i].index=i;
          samples[i].val=distribution(generator);
        }

        //cout<<"Samples created"<<endl;

        sort(samples.begin(),samples.end(),markedSampComp);
        //cout<<"Samples sorted"<<endl;

        vector<int> pickedCoords (M,0);
        float runningProbSum = 0.0;
        int sampInd=0,coordInd=0;

        while((sampInd<M)&&(coordInd<numOfLocs1+numOfLocs2)){
          while((runningProbSum<samples[sampInd].val)&&(coordInd<numOfLocs1+numOfLocs2)){
            runningProbSum += probs[coordInd]/probSum;
            coordInd++;
          }
          while((samples[sampInd].val<=runningProbSum)&&(sampInd<M)&&(coordInd<numOfLocs1+numOfLocs2)){
            pickedCoords[samples[sampInd].index]=coordInd;
            sampInd++;
          }
        }
        //cout<<"Locations picked"<<endl;

        for(int i=0;(i<M)&&(total<(rank1+rank2));i++)
          if((point[pickedCoords[i]]<1)&&(chosen[pickedCoords[i]]==false)){
            //cout<<pickedCoords[i]<<" "<<score(pickedCoords[i])<<" "<<exp(eps*score(pickedCoords[i])/2)<<" "<<tolerance<<" "<<((tolerance*probs[i])+1)<<" "<<probs[i]<<endl;
            if (distribution(generator) <= (exp(eps*score(pickedCoords[i])/2)/((probs[i])+0.01))){
              point[pickedCoords[i]] += eta;
              chosen[pickedCoords[i]]=true;
              total++;
              //cout<<"BUMP"<<" "<<score(pickedCoords[i])<<" "<<probs[i]<<" "<<total<<endl;
            }
          }
      }
    }
}

bool continuousGreedy::markedSampComp(markedSamp s1, markedSamp s2){
  if (s1.val<s2.val)
    return true;
  else
    return false;
}

float continuousGreedy::algoWrapper(){
  exponentialMechanism();
  ctsScore = 0;

  for(int i = 0; i<numOfSamples; i++){
    set<int, less<int> > sampleSet;
    for(int j = 0; j<(numOfLocs1 + numOfLocs2);j++)
      if (hashes[j][i]<point[j])
        sampleSet.insert(j);
    ctsScore += funValue(sampleSet);
  }

  ctsScore /= numOfSamples;

  returnSet = swapRounding();

  return (funValue(returnSet));

}

std::set<int, less<int> > continuousGreedy::swapRounding(){

  std::set<int, less<int> > ans, part2;
  std::set<int, less<int> >::iterator itr;
  ans = partitionRound(point, numOfLocs1);
  float* part2Begin = &point[numOfLocs1];
  if (numOfLocs2>0){
    part2 = partitionRound(part2Begin, numOfLocs2);
    for(itr = part2.begin(); itr!= part2.end(); itr++)
      ans.insert((*itr)+numOfLocs1);
    }
  return ans;

}

std::set<int, less<int> > continuousGreedy::partitionRound(float* point, int n){
  //Do not touch i and j.
  int i = 0, iCheck = 0, j = 0, jCheck = 0;
  float checkSum = 0;
  float lowTolerance = 0.01, hiTolerance=0.99;
  std::vector<float> probs(2);
  float sum = 0;
  int winner;

  i = 0;
  std::set<int, less<int> > final;
  do{
    //Rounding for clarity
    for(int k = 0; k<n; k++)
      if (point[k]<lowTolerance)
        point[k] = 0;
      else if (point[k]>hiTolerance)
        point[k] = 1;

    do{
      while ((i<n)&&(point[i]<lowTolerance))
        i++;
      iCheck = i;
      while ((i<n)&&(point[i]>hiTolerance)){
        final.insert(i);
        point[i] = 0;
        i++;
      }
    }while((i!=iCheck)&&(i<n));

    if (i==n)
      break;
    //At this point 0<point[i]<1
    j = i+1;
    do{
      while ((j<n)&&(point[j]<lowTolerance))
        j++;
      jCheck = j;
      while ((j<n)&&(point[j]>hiTolerance)){
        point[j] = 0;
        final.insert(j);
        j++;
      }
    }while((j!=jCheck)&&(j<n));

    if (j==n)
      break;
    //At this point 0<point[j]<1

    sum = point[i] + point[j];
    if (sum>1){
      probs[0] = (1-point[j])/(2-sum);
      probs[1] = (1-point[i])/(2-sum);
    }
    else{
      probs[0] = point[i]/sum;
      probs[1] = point[j]/sum;
    }
    std::discrete_distribution<> distributionScore(probs.begin(), probs.end());
    winner = distributionScore(generator);
    if (sum>1)
      if (!winner){
        final.insert(i);
        point[i] = 0;
        point[j] = sum - 1;
        i = j;
      }
      else{
        final.insert(j);
        point[j] = 0;
        point[i] = sum - 1;
      }
    else{
      if (!winner){
        point[i] = sum;
        point[j] = 0;
      }
      else{
        point[i] = 0;
        point[j] = sum;
        i = j;
      }
    }
}while(i<n);

return final;

}

set<int, less <int> > continuousGreedy::exponentialMechanism2(){
  int check1=0 , check2=0;

  scoreComputer();

  set<int, less <int> > basisSet;
  set<int, less <int> >::iterator itr;
  int chosen;

  while((check1+check2)<(rank1+rank2)){
    //Initialise the probabilities.
    std::vector<float> probs;
    for(int i=0;i<(numOfLocs1 + numOfLocs2);i++)
      probs.push_back(exp(eps*(scores[i])/2));
    for(auto itr: basisSet)
      probs[itr]=0;
    //Remove elements that make the set lie outside the matroid.
    if (check1==rank1)
      for(int i=0;i<numOfLocs1;i++)
        probs[i] = 0;
    if (check2==rank2)
      for(int i=numOfLocs1;i<(numOfLocs1+numOfLocs2);i++)
        probs[i] = 0;
    //Choose an element and update the matroid budgets.
    std::discrete_distribution<> distributionScore(probs.begin(), probs.end());
    chosen = distributionScore(generator);
    if (chosen<numOfLocs1)
      check1++;
    else
      check2++;
    basisSet.insert(chosen);
    point[chosen]+=(1.0/T);

  }

  return basisSet;
}

float continuousGreedy::algoWrapper2(){
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

  std::vector <set<int, less <int> > > chosenSets(T);
  set<int, less <int> > set1;

  for(int i=0;i<(numOfLocs1 + numOfLocs2);i++)
    set1.insert(i);

  set<int, less <int> >::iterator itr;
  int j=0;

  for(int i=0;i<T;i++){
    chosenSets[i] = exponentialMechanism2();

    float tempSum = 0;
    for(int j=0;j<(numOfLocs1 + numOfLocs2);j++)
      tempSum += point[j];
    if ((tempSum> ((rank1+rank2)*float(i+1)/T +0.02))||(tempSum< (rank1+rank2)*float(i+1)/T-0.02))
      cout<<"PROBLEM"<<tempSum<<' '<<(rank1+rank2)*float(i+1)/T<<endl;

  }

  ctsScore = 0;

  if (!(direct)){
  for(int i = 0; i<numOfSamples; i++){
    set<int, less<int> > sampleSet;
    for(int j = 0; j<(numOfLocs1 + numOfLocs2);j++)
      if (hashes[j][i]<point[j])
        sampleSet.insert(j);
    ctsScore += funValue(sampleSet);
  }
  }

  ctsScore /= numOfSamples;

  returnSet = swapRounding();

  return (funValue(returnSet));

}

#endif  // DPO_CG_H
