#ifndef DPO_INITIALISE_H
#define DPO_INITIALISE_H

#include "preproc.hpp"

#include<iostream>
#include<fstream>
#include<functional>
#include<random>
#include<set>
#include<cmath>
#include<chrono>
#include<string>

using namespace std;

class initialise{
  public:
    std::default_random_engine generator;
    vector <float> locLat, locLon; //Location coordinates
    float *lat, *lon; //Data coordinates
    int numOfLocs1, numOfLocs2, numOfData;
    float latNorth, latSouth, latEast, latWest;
    float lonNorth, lonSouth, lonEast, lonWest;
    vector <vector <int>> orderedLocations; //For each datapoint stores
                                              //locations in increasing order of
                                              //distance from that point
    vector<int> increasingOrder;              //Scratchpad for noting locations
                                              //in increasing order of distance
                                              //from a fixed datapoint.
    vector <float> distances;                 //distances[i][j] stores distance
                                              //of ith point with jth location
    vector<vector <float>> allDistances;

    float rank1, rank2;
    float eps, delta;
    std::set <int, std::less <int> > final;
    float e = 2.7182818284;
    float norm;

    ~initialise();

    initialise  (int givenNumOfLocs1,
      int givenNumOfLocs2,
      int givenNumOfData,
      float givenEps,
      float givenDelta,
      int givenRank1,
      int givenRank2,
      vector<float> box,
      vector<float> givenLocLat,
      vector<float> givenLocLon);

    void drawData();
    void hashRenew();

    float funValue(std::set<int, std::less <int> > argSet);

    void sortLocations(int,int);

    void setPrint(std::set<int, std::less <int> > argSet);

    unsigned long long rdtsc();

};

initialise::~initialise(){
  //Release all dynamically allocated memory.
  delete[] lat;
  delete[] lon;
}

initialise::initialise  (int givenNumOfLocs1,
  int givenNumOfLocs2,
  int givenNumOfData,
  float givenEps,
  float givenDelta,
  int givenRank1,
  int givenRank2,
  vector<float> box,
  vector<float> givenLocLat,
  vector<float> givenLocLon){
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

  /*
  Variables initialised in order.
  1. std::default_random_engine generator;
  2. int numOfLocs1, int numOfLocs2, numOfData;
  3. float *locLat, *locLon; //Location coordinates
  4. float *lat, *lon; //Data coordinates
  5. float rank1, rank2;
  6. int T;
  7. float eta;
  8. float **hashes;
  9. float eps;
  10. float *scores;
  11. float *points;
  12. float norm;

  */

  //1
  generator.seed(rdtsc());
  //2
  numOfLocs1 = givenNumOfLocs1, numOfLocs2 = givenNumOfLocs2,
  numOfData = givenNumOfData;

  //3

  float latNorth = box[0];
  float latSouth = box[1];
  float latEast = box[2];
  float latWest = box[3];
  float lonNorth = box[4];
  float lonSouth = box[5];
  float lonEast = box[6];
  float lonWest = box[7];

  //12
  float max = (fabs(lonNorth - lonSouth) + fabs(latNorth-latSouth));
  if (max<(fabs(lonEast - lonWest) + fabs(latEast - latWest)))
    max = (fabs(lonEast - lonWest) + fabs(latEast - latWest));
  if (max<(fabs(lonNorth-lonWest) + fabs(latNorth - latWest)))
    max = (fabs(lonNorth-lonWest) + fabs(latNorth - latWest));
  if (max < (fabs(lonEast-lonSouth) + fabs(latEast - latSouth)))
    max = (fabs(lonEast-lonSouth) + fabs(latEast - latSouth));
  norm = max;

  locLat = givenLocLat;
  locLon = givenLocLon;

  //4
  fstream f;
  string line;
  lat = new float[numOfData];
  lon = new float[numOfData];

  distances.resize(numOfLocs1 + numOfLocs2,0);
  for(int i=0;i<numOfLocs1+numOfLocs2; i++)
    increasingOrder.push_back(i);

  f.open("../data/uber-data-small.csv");
  for(int i = 0; i<numOfData; i++){
    getline(f,line,',');
    lat[i] = stof(line);
    getline(f,line);
    lon[i] = stof(line);
    allDistances.push_back(distances);
    float min = norm;
    for(int j=0;j<numOfLocs1+numOfLocs2;j++){
      float dist = (fabs(lat[i] - locLat[j]) + fabs(lon[i] - locLon[j]))/norm;
      distances[j] = dist;
      allDistances[i][j] = dist;
      //cout<<"Distances for point"<<lat[i]<<","<<lon[i]<<" location "<<locLat[j]<<","<<locLon[j]<<" is "<<distances[j]<<endl;
    }

    sortLocations(0,numOfLocs1+numOfLocs2-1);
    orderedLocations.push_back(increasingOrder);
    //for(int i=0;i<numOfLocs1+numOfLocs2;i++)
    //  cout<<increasingOrder[i]<<";"<<distances[increasingOrder[i]]<<"\t";
  }
  f.close();

  //5
  rank1 = givenRank1;
  rank2 = givenRank2;

  //9
  eps = givenEps;
  delta = givenDelta;
}

void initialise::drawData(){
  //Do not add dataInit here, data shuffled uniformly before this point
  fstream f;
  string line;
  orderedLocations.clear();

  f.open("../data/uber-data-small.csv");
  for(int i = 0; i<numOfData; i++){
    getline(f,line,',');
    lat[i] = stof(line);
    getline(f,line);
    lon[i] = stof(line);

    float min = norm;
    for(int j=0;j<numOfLocs1+numOfLocs2;j++){
      float dist = (fabs(lat[i] - locLat[j]) + fabs(lon[i] - locLon[j]))/norm;
      distances[j] = dist;
      allDistances[i][j] = dist;
      //cout<<"Distances for point"<<lat[i]<<","<<lon[i]<<" location "<<locLat[j]<<","<<locLon[j]<<" is "<<distances[j]<<endl;
    }

    sortLocations(0,numOfLocs1+numOfLocs2-1);
    orderedLocations.push_back(increasingOrder);
    //for(int i=0;i<numOfLocs1+numOfLocs2;i++)
    //  cout<<increasingOrder[i]<<";"<<distances[increasingOrder[i]]<<"\t";
  }
  f.close();
}


void initialise::sortLocations(int start,int end){
  int len = end-start + 1;
  if(len>1){
    sortLocations(start,start+(len/2)-1);
    sortLocations(start+(len/2),end);
    vector<int> temp(len,0);

    int i=start, j=start+len/2, t=0;

    for(;i<start+(len/2)&&(j<=end);)
      if (distances[increasingOrder[i]]>distances[increasingOrder[j]])
        temp[t++] = increasingOrder[j++];
      else
        temp[t++] = increasingOrder[i++];

    for(;i<start+(len/2);)
      temp[t++] = increasingOrder[i++];

    for(;j<=end;)
      temp[t++] = increasingOrder[j++];

    for(t=0;t<len;t++)
      increasingOrder[start+t] = temp[t];
  }
}

float initialise::funValue(set<int, less <int> > argSet){
  /*Computes:
  numOfData - \sum_{i \in Data} min_{j \in Locations} (Manhattan-Dist(i,j)/0.266)
  locLat and locLon are coordinates of Locations, lat and lon are coordinates
  of Data.
  */
  if (!(argSet.size()))
    return (0);
  else{
    float minLat[numOfData], minLon[numOfData];
    float min = 1;
    float sum = 0;
    for(int j = 0;j<numOfData; j++){
      min = -1;
      for(auto itr:argSet){
        float val = (fabs(locLat[itr] - lat[j]) + fabs(locLon[itr] - lon[j]))/(norm);
        //val is Manhattan Distance between jth data point and *itr-th location
        if (min<0){
          min = val;
          minLat[j] = locLat[itr];
          minLon[j] = locLon[itr];
        }
        else if (val<min){
          min = val;
          minLat[j] = locLat[itr];
          minLon[j] = locLon[itr];
        }
      }
      //min is min_{*itr \in Locations} Manhattan Distance(*itr,j)
      if(min>=0)
        sum+= min;
    }
    //sum is \sum_{j \in Data} min_{*itr \in Locations} Manhattan Distance(*itr,j)

    return numOfData - sum;
  }
}

void initialise::setPrint(set<int, less <int> > argSet){
  set<int, less <int> >::iterator i;
  cout<<"Set: ";
  for(auto i:argSet)
    cout<<i<<' ';
}

unsigned long long initialise::rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}

#endif  // DPO_CG_H
