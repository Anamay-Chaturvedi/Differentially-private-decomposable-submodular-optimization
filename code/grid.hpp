#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

const float pi = 3.1415;

float twist1(float input){
  return pow(sin((pi/2)*input),1);
}

float twist2(float input){
  return pow(sin((pi/2)*input),2);
}

float twister(int which, int howMany, float what){
  if (which==1)
    for(int i=0;i<howMany;i++)
      what = twist1(what);
  else if (which==2)
    for(int i=0;i<howMany;i++)
      what = twist2(what);
  return what;
}

void createGrid(int numOfLocs, int gridCols, vector <float> box, int spurious){
  //float latNorth=40.8, latSouth=40.743775, latEast=40.79, latWest=40.756412;
  //float lonNorth=-73.97, lonSouth=-73.973717, lonEast=-73.94, lonWest=-74.003862;
  //float latNorth = 40.814381, latSouth = 40.743775, latEast = 40.743775, latWest = 40.814381;
  //float lonNorth = -73.931548, lonSouth = -74.003862, lonEast = -73.931548, lonWest = -74.003862;
  //float latNorth = 40.81794, latWest = 40.71315, latSouth = 40.7, latEast = 40.80204;
  //float lonNorth = -73.96483, lonWest = -74.03, lonSouth = -73.99197, lonEast = -73.91436;
  float latNorth = box[0];
  float latSouth = box[1];
  float latEast = box[2];
  float latWest = box[3];
  float lonNorth = box[4];
  float lonSouth = box[5];
  float lonEast = box[6];
  float lonWest = box[7];

  float locLat[numOfLocs], locLon[numOfLocs];

  /*
  for(int i = 0;i<numOfLocs-spurious;i++){
    locLat[i] = latSouth + (latWest - latSouth)*twister(2,1,(i/gridCols)*(1.0/((numOfLocs-1)/gridCols)))
                  + (latEast - latSouth)*twister(1,4,((i%gridCols)*(1.0/(gridCols-1.0))));
    locLon[i] = lonSouth + (lonWest - lonSouth)*twister(2,1,(i/gridCols)*(1.0/((numOfLocs-1)/gridCols)))
                  + (lonEast - lonSouth)*twister(1,4,((i%gridCols)*(1.0/(gridCols-1.0))));
  }
  */

  numOfLocs = numOfLocs - spurious;

  for(int i = 0;i<numOfLocs;i++){
    locLat[i] = latSouth + (latWest - latSouth)*(i/gridCols)*(1.0/((numOfLocs-1)/gridCols))
                  + (latEast - latSouth)*((i%gridCols)*(1.0/(gridCols-1.0)));
    locLon[i] = lonSouth + (lonWest - lonSouth)*(i/gridCols)*(1.0/((numOfLocs-1)/gridCols))
                  + (lonEast - lonSouth)*((i%gridCols)*(1.0/(gridCols-1.0)));
  }
  for(int i = numOfLocs;i<numOfLocs+spurious;i++){
    locLat[i] = latNorth;
    locLon[i] = lonNorth;
  }

  fstream f;
  f.open("twistGrid.csv",ios::out);
  for(int i=0;i<numOfLocs+spurious;i++)
    f<<locLat[i]<<','<<locLon[i]<<endl;
  f.close();


}
