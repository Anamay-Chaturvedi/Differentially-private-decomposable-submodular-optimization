#include<iostream>
#include<string>
#include<vector>
#include<set>
#include<cmath>
#include<fstream>
#include<chrono>
#include<ctime>

#include "algorithms/greedy.hpp"
#include "algorithms/continuousGreedy.hpp"
#include "algorithms/diffPrivateGreedy.hpp"
#include "algorithms/randomSet.hpp"
#include "preproc.hpp"


using namespace std;

void experiment1(int);                  //The vector<int> rank is the set of ranks the experiment is run over.
void experiment2(int);                  //The vector<int> dataSizes is the set of dataSizes the experiment is run over.
void experiment3(int);

void experimentTemplate(int numOfLocs1, //Number of locations in the first partition
  int numOfLocs2,                       //Number of locations in the second partition
  int numberOfDataPoints,               //Number of pick-ups used to define the location utility
  int numberOfSamples,                  //Number of random vectors used by the continuous greedy to construct the multilinear approximation
  int T,                                //Number of rounds used in the continuous greedy
  bool guptaPrivacy,                    //If true the Gupta privacy parameter is used, if false, the best composition law is used
  float epsilon,                        //The privacy parameter all algorithms are required to fulfill.
  int rank1,                            //Number of locations that may be picked from the first matroid.
  int rank2,                            //Number of locations that may be picked from the second matroid.
  int numOfRuns,                        //Number of times the experiment is conducted. Data is drawn afresh every 10 runs.
  vector<float> &means,                 //A vector used to record the empirical mean.
  vector<float> &stdDevs,               //A vector used to record the empirical standard deviation.
  vector <float> box,                   //A vector used to pass the coordinates of the bounding box in Manhattan.
  vector <float> locLat,                //A vector used to pass the coordinates of the grid of locations to choose from.
  vector <float> locLon,                //A vector used to pass the coordinates of the grid of locations to choose from.
  bool direct,
  bool partition);

int main(){
  int numOfRuns = 400;
  //cout<<endl<<"Enter number of runs for the cardinality constraint experiment:"<<endl;
  //cin>>numOfRuns;
  auto start1 = std::chrono::system_clock::now();
  experiment1(numOfRuns);
  auto end1 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds1 = end1-start1;
  cout<<"Experiment 1 "<<elapsed_seconds1.count()<<endl;

  numOfRuns = 400;
  //cout<<endl<<"Enter number of runs for the cardinality constraint experiment:"<<endl;
  //cin>>numOfRuns;
  auto start2 = std::chrono::system_clock::now();
  experiment2(numOfRuns);
  auto end2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds2 = end2-start2;
  cout<<"Experiment 2 "<<elapsed_seconds2.count()<<endl;

  return 0;
}

void experimentTemplate(int numOfLocs1,
  int numOfLocs2,
  int numberOfDataPoints,
  int numberOfSamples,
  int T,
  bool guptaPrivacy,
  float epsilon,
  int rank1,
  int rank2,
  int numOfRuns,
  vector<float> &means,
  vector<float> &stdDevs,
  vector <float> box,
  vector <float> locLat,
  vector <float> locLon,
  bool direct,
  bool partition){

  const float e = 2.7182818284;
  float delta = 1.0/pow(numberOfDataPoints,1.5);
  cout<<"continuousGreedy's eps is "<<2*log(1 + epsilon/(2*(3 - log(delta))))<<" and delta is "<<delta<<endl;

  dataInit(box);

  std::vector<float>  greedyUtils(numOfRuns,0), expMechUtils(numOfRuns,0),
                      randomUtils(numOfRuns,0), continuousGreedyUtils(numOfRuns,0);

  float greedyUtilsAvg = 0, greedyUtilsStdDev = 0, expMechUtilsAvg = 0,
        expMechUtilsStdDev = 0, randomUtilsAvg = 0, randomUtilsStdDev = 0,
        continuousGreedyUtilsAvg = 0, continuousGreedyUtilsStdDev = 0,
        ctsAvg = 0, ctsStdDev = 0;

  cout<<"numOfLocs1"<<' '<<"numOfLocs2"<<' '<<"numberOfDataPoints"<<' '<<"numberOfSamples"<<' '<<
                    "T"<<' '<<"epsilon"<<' '<<"rank1"<<' '<<"rank2"<<endl;
  cout<<numOfLocs1<<' '<<numOfLocs2<<' '<<numberOfDataPoints<<' '<<numberOfSamples<<' '<<
                    T<<' '<<epsilon<<' '<<rank1<<' '<<rank2<<endl;

  continuousGreedy A1 (numOfLocs1, numOfLocs2, numberOfDataPoints, numberOfSamples,
                    T, epsilon, delta, rank1, rank2, box, locLat, locLon, direct);
  greedy A2 (numOfLocs1, numOfLocs2, numberOfDataPoints,
                    epsilon, 0, rank1, rank2, box, locLat, locLon);
  diffPrivateGreedy A3 (numOfLocs1, numOfLocs2, numberOfDataPoints,
                    epsilon, delta, rank1, rank2, box, guptaPrivacy, locLat, locLon);

  randomSet A4 (numOfLocs1, numOfLocs2, numberOfDataPoints,
                    epsilon, 0, rank1, rank2, box, locLat, locLon);

  for(int expIndex  = 0; expIndex<numOfRuns; expIndex++){
    if (!(expIndex%10)){
      dataInit(box); //Shuffles and writes into Uber-data-small
      A1.drawData();
      A2.drawData();
      A3.drawData();
      A4.drawData();

      cout<<"Data redrawn"<<endl;
      cout<<"i\t CG\t Cts\t Gdy\t R\t DPG\t CGavg \t DPGavg\n";

    }
    if (direct==false)
      A1.hashRenew();

    if (!(partition)){
      continuousGreedyUtils[expIndex] = A1.algoWrapper();
    }
    else{
      continuousGreedyUtils[expIndex] = A1.algoWrapper2();
    }
    continuousGreedyUtilsAvg += continuousGreedyUtils[expIndex];
    continuousGreedyUtilsStdDev += continuousGreedyUtils[expIndex]*continuousGreedyUtils[expIndex];

    greedyUtils[expIndex] = A2.algoWrapper();
    greedyUtilsAvg += greedyUtils[expIndex];
    greedyUtilsStdDev += greedyUtils[expIndex]*greedyUtils[expIndex];

    expMechUtils[expIndex] = A3.algoWrapper();
    expMechUtilsAvg += expMechUtils[expIndex];
    expMechUtilsStdDev += expMechUtils[expIndex]*expMechUtils[expIndex];

    randomUtils[expIndex] = A4.algoWrapper();
    randomUtilsAvg += randomUtils[expIndex];
    randomUtilsStdDev += randomUtils[expIndex]*randomUtils[expIndex];


    cout<<expIndex<<"\t"<<continuousGreedyUtils[expIndex]<<"\t"<<A1.ctsScore<<"\t"<<greedyUtils[expIndex]
        <<"\t"<<randomUtils[expIndex]<<"\t"<<expMechUtils[expIndex]<<"\t"<<continuousGreedyUtilsAvg/(expIndex+1)
        <<"\t"<<expMechUtilsAvg/(expIndex+1)<<endl;

  }

  continuousGreedyUtilsAvg /= numOfRuns;
  continuousGreedyUtilsStdDev /= numOfRuns;
  greedyUtilsAvg /= numOfRuns;
  greedyUtilsStdDev /= numOfRuns;
  expMechUtilsAvg/= numOfRuns;
  expMechUtilsStdDev/= numOfRuns;
  randomUtilsAvg/= numOfRuns;
  randomUtilsStdDev/= numOfRuns;
  cout<<endl<<greedyUtilsStdDev<<endl;
  cout<<endl<<greedyUtilsAvg*greedyUtilsAvg<<endl;

  means[0] = continuousGreedyUtilsAvg;
  stdDevs[0] = sqrt(continuousGreedyUtilsStdDev - continuousGreedyUtilsAvg*continuousGreedyUtilsAvg);
  means[1] = greedyUtilsAvg;
  stdDevs[1] = sqrt(greedyUtilsStdDev - greedyUtilsAvg*greedyUtilsAvg);
  means[2] = expMechUtilsAvg;
  stdDevs[2] = sqrt(expMechUtilsStdDev - expMechUtilsAvg*expMechUtilsAvg);
  means[3] = randomUtilsAvg;
  stdDevs[3] = sqrt(randomUtilsStdDev - randomUtilsAvg*randomUtilsAvg);

  cout<<"Average continouous greedy Utility: "<<means[0]<<" and standard deviation: "<<stdDevs[0]<<endl;
  cout<<"Average greedy Utility: "<<means[1]<<" and standard deviation: "<<stdDevs[1]<<endl;
  cout<<"Average discrete private greedy Utility: "<<means[2]<<" and standard deviation: "<<stdDevs[2]<<endl;
  cout<<"Average random Utility: "<<means[3]<<" and standard deviation: "<<stdDevs[3]<<endl;
}

void experiment1(int numOfRuns){
  float latNorth = 40.82, latWest = 40.71315, latSouth = 40.7, latEast = 40.80;
  float lonNorth = -73.96, lonWest = -74.03, lonSouth = -73.99197, lonEast = -73.93;

  vector<float> box = {latNorth,latSouth,latEast,latWest,lonNorth,lonSouth,lonEast,lonWest};

  vector <float> ranks = {10,13,16,19,22,25};  //This is the set of ranks that the experiment will be run on.
  int numOfExp = ranks.size();
  int numOfLocs = 100;
  int numOfSpurious = 80;
  int gridCols = 5;

  vector<vector<float> > means;
  vector<vector<float> > stdDevs;

  vector<float> locLat(numOfLocs);
  vector<float> locLon(numOfLocs);


  for(int i = 0;i<numOfLocs-numOfSpurious;i++){
    locLat[i] = latSouth + (latWest - latSouth)*(i/gridCols)*(1.0/((numOfLocs-1)/gridCols))
                  + (latEast - latSouth)*((i%gridCols)*(1.0/(gridCols-1.0)));
    locLon[i] = lonSouth + (lonWest - lonSouth)*(i/gridCols)*(1.0/((numOfLocs-1)/gridCols))
                  + (lonEast - lonSouth)*((i%gridCols)*(1.0/(gridCols-1.0)));
  }
  for(int i = numOfLocs-numOfSpurious;i<numOfLocs;i++){
    locLat[i] = latNorth;
    locLon[i] = lonNorth;
  }

  for(int i=0;i<numOfExp;i++){
    std::vector<float> tempVec(4);
    means.push_back(tempVec);
    stdDevs.push_back(tempVec);
    experimentTemplate(numOfLocs,//int numOfLocs1,
      0,//int numOfLocs2,
      100,//int numberOfDataPoints,
      1000,//int numberOfSamples,
      3,//int T,
      false,
      0.1,//float epsilon,
      ranks[i],//int rank1,
      0,//int rank2,
      numOfRuns,//int numOfRuns
      means[i],
      stdDevs[i],
      box,
      locLat,
      locLon,
      true,
      false);
  }

  cout<<"Ranks"<<endl;
  for(int i=0;i<numOfExp;i++)
    cout<<ranks[i]<<',';
  cout<<endl;
  cout<<"Continuous greedy"<<endl;
  cout<<"Means: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<means[i][0]<<',';
  cout<<"]"<<endl;
  cout<<"Standard deviations: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<stdDevs[i][0]<<',';
  cout<<"]"<<endl;

  cout<<"Greedy: "<<endl;

  cout<<"Means: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<means[i][1]<<',';
  cout<<"]"<<endl;
  cout<<"Standard deviations: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<stdDevs[i][1]<<',';
  cout<<"]"<<endl;


  cout<<"Discrete private greedy"<<endl;

  cout<<"Means: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<means[i][2]<<',';
  cout<<"]"<<endl;
  cout<<"Standard deviations: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<stdDevs[i][2]<<',';
  cout<<"]"<<endl;


  cout<<"Random"<<endl;

  cout<<"Means: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<means[i][3]<<',';
  cout<<"]"<<endl;
  cout<<"Standard deviations: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<stdDevs[i][3]<<',';
  cout<<"]"<<endl;

  fstream file1("../results/exp1.csv", ios::out);
  file1<<"Ranks,CGavg,CGstdDev,Gavg,GstdDev,DPGavg,DPGstdDev,Ravg,RstdDev,"<<numOfRuns<<endl;
  for(int i=0;i<numOfExp;i++){
    file1<<ranks[i];
    for(int j=0;j<4;j++)
      file1<<','<<means[i][j]<<','<<stdDevs[i][j];
    file1<<"\n";
  }
  file1.close();
}

void experiment2(int numOfRuns){
  float latNorth=40.814381, latSouth=40.743775, latEast=40.743775, latWest=40.814381;
  float lonNorth=-73.931548, lonSouth=-74.003862, lonEast=-73.931548, lonWest=-74.003862;

  vector<float> box = {latNorth,latSouth,latEast,latWest,lonNorth,lonSouth,lonEast,lonWest};

  vector <float> dataSizes = {1000,2000,3000,4000,5000,6000,7000,8000,9000,10000};  //This is the set of datasizes that the experiment will be run on.
  int numOfExp = dataSizes.size();
  cout<<numOfExp;

  vector<vector<float> > means(numOfExp);
  vector<vector<float> > stdDevs(numOfExp);

  float lat0 = 40.75, lon0 = -74.0025;
  float lat1 = 40.7525, lon1 = -73.9975;
  float lat2 = 40.7825, lon2 = -73.9575;
  //1 should be best, 0 should be similar but worse, and 2 should provide the most marginal benefit
  vector <float> locLat={lat0, lat1, lat2};
  vector <float> locLon={lon0, lon1, lon2};

  for(int i=0;i<numOfExp;i++){
    vector<float> tempVec(10);
    means[i] = tempVec;
    vector<float> tempVec2(10);
    stdDevs[i] = tempVec2;
    experimentTemplate(1,//int numOfLocs1,
      2,//int numOfLocs2,
      dataSizes[i],//int numberOfDataPoints,
      1000,//int numberOfSamples,
      7,//int T,
      true,
      0.1,//float epsilon,
      1,//int rank1,
      1,//int rank2,
      numOfRuns,//int numOfRuns
      means[i],
      stdDevs[i],
      box,
      locLat,
      locLon,
      true,
      true);
  }

  cout<<"DataSizes"<<endl;
  for(int i=0;i<numOfExp;i++)
    cout<<dataSizes[i]<<',';
  cout<<endl;
  cout<<"Continuous greedy"<<endl;
  cout<<"Means: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<means[i][0]<<',';
  cout<<"]"<<endl;
  cout<<"Standard deviations: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<stdDevs[i][0]<<',';
  cout<<"]"<<endl;

  cout<<"Greedy: "<<endl;

  cout<<"Means: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<means[i][1]<<',';
  cout<<"]"<<endl;
  cout<<"Standard deviations: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<stdDevs[i][1]<<',';
  cout<<"]"<<endl;


  cout<<"Discrete private greedy"<<endl;

  cout<<"Means: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<means[i][2]<<',';
  cout<<"]"<<endl;
  cout<<"Standard deviations: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<stdDevs[i][2]<<',';
  cout<<"]"<<endl;


  cout<<"Random"<<endl;

  cout<<"Means: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<means[i][3]<<',';
  cout<<"]"<<endl;
  cout<<"Standard deviations: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<stdDevs[i][3]<<',';
  cout<<"]"<<endl;

  fstream file2("../results/exp2.csv", ios::out);
  file2<<"DataSizes,CGavg,CGstdDev,Gavg,GstdDev,DPGavg,DPGstdDev,Ravg,RstdDev,"<<numOfRuns<<endl;
  for(int i=0;i<numOfExp;i++){
    file2<<dataSizes[i];
    for(int j=0;j<4;j++)
      file2<<','<<means[i][j]<<','<<stdDevs[i][j];
    file2<<"\n";
  }
  file2.close();

}

void experiment3(int numOfRuns){
  float latNorth=40.814381, latSouth=40.743775, latEast=40.743775, latWest=40.814381;
  float lonNorth=-73.931548, lonSouth=-74.003862, lonEast=-73.931548, lonWest=-74.003862;

  vector<float> box = {latNorth,latSouth,latEast,latWest,lonNorth,lonSouth,lonEast,lonWest};

  vector <float> epsilons = {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45};  //This is the set of datasizes that the experiment will be run on.
  int numOfExp = epsilons.size();
  cout<<numOfExp;

  vector<vector<float> > means(numOfExp);
  vector<vector<float> > stdDevs(numOfExp);

  float lat0 = 40.75, lon0 = -74.0025;
  float lat1 = 40.7525, lon1 = -73.9975;
  float lat2 = 40.7825, lon2 = -73.9575;
  //1 should be best, 0 should be similar but worse, and 2 should provide the most marginal benefit
  vector <float> locLat={lat0, lat1, lat2};
  vector <float> locLon={lon0, lon1, lon2};

  for(int i=0;i<numOfExp;i++){
    vector<float> tempVec(10);
    means[i] = tempVec;
    vector<float> tempVec2(10);
    stdDevs[i] = tempVec2;
    experimentTemplate(1,//int numOfLocs1,
      2,//int numOfLocs2,
      5000,//int numberOfDataPoints,
      1000,//int numberOfSamples,
      7,//int T,
      true,
      epsilons[i],//float epsilon,
      1,//int rank1,
      1,//int rank2,
      numOfRuns,//int numOfRuns
      means[i],
      stdDevs[i],
      box,
      locLat,
      locLon,
      true,
      true);
  }

  cout<<"Epsilons"<<endl;
  for(int i=0;i<numOfExp;i++)
    cout<<epsilons[i]<<',';
  cout<<endl;
  cout<<"Continuous greedy"<<endl;
  cout<<"Means: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<means[i][0]<<',';
  cout<<"]"<<endl;
  cout<<"Standard deviations: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<stdDevs[i][0]<<',';
  cout<<"]"<<endl;

  cout<<"Greedy: "<<endl;

  cout<<"Means: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<means[i][1]<<',';
  cout<<"]"<<endl;
  cout<<"Standard deviations: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<stdDevs[i][1]<<',';
  cout<<"]"<<endl;


  cout<<"Discrete private greedy"<<endl;

  cout<<"Means: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<means[i][2]<<',';
  cout<<"]"<<endl;
  cout<<"Standard deviations: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<stdDevs[i][2]<<',';
  cout<<"]"<<endl;


  cout<<"Random"<<endl;

  cout<<"Means: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<means[i][3]<<',';
  cout<<"]"<<endl;
  cout<<"Standard deviations: "<<endl<<"[";
  for(int i=0;i<numOfExp;i++)
    cout<<stdDevs[i][3]<<',';
  cout<<"]"<<endl;

  fstream file2("../results/exp3.csv", ios::out);
  file2<<"DataSizes,CGavg,CGstdDev,Gavg,GstdDev,DPGavg,DPGstdDev,Ravg,RstdDev,"<<numOfRuns<<endl;
  for(int i=0;i<numOfExp;i++){
    file2<<epsilons[i];
    for(int j=0;j<4;j++)
      file2<<','<<means[i][j]<<','<<stdDevs[i][j];
    file2<<"\n";
  }
  file2.close();

}
