#ifndef DPO_PREPROC_H
#define DPO_PREPROC_H

#include<iostream>
#include<fstream>
#include<random>
#include<vector>
#include<string>
#include<algorithm>
#include<chrono>

using namespace std;

void dataInit(vector <float> box){

  float latNorth = box[0];
  float latSouth = box[1];
  float latEast = box[2];
  float latWest = box[3];
  float lonNorth = box[4];
  float lonSouth = box[5];
  float lonEast = box[6];
  float lonWest = box[7];

  long num = 0;
  float prob = 0;
  long numOfData = 10000;
  int index;

  fstream f1 ("../data/uber-raw-data-apr14.csv"), f2 ("../data/uber-data-small.csv", ios::out);

  string line;
  vector <string> file(numOfData);

  float lat, lon;

  std::vector<int> indices;

  for(int i=0;i<numOfData;i++)
    indices.push_back(1);

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution1(0.0,1.0);
  std::discrete_distribution<> distribution2(indices.begin(), indices.end());
  generator.seed(time(NULL));

  //Create a grid of locations with gridCols many columns.
  //float latNorth = 40.81794, latWest = 40.71315, latSouth = 40.6866, latEast = 40.80204;
  //float lonNorth = -73.96483, lonWest = -74.04519, lonSouth = -73.99197, lonEast = -73.91436;

  //North 40.81794, -73.96483 done
  //West 40.71315, -74.04519 done
  //East 40.80204, -73.91436 done
  //South 40.6866, -73.99197 done

  getline(f1,line);

  while(!(f1.eof())){
    getline(f1,line,',');
    if (line=="")
      break;
    getline(f1,line,',');
    //cout<<" "<<line<<", ";
    lat = stof(line);
    getline(f1,line,',');
    //cout<<line<<endl;
    lon = stof(line);
    getline(f1,line);

    int check1 = 0, check2 = 0, check3 = 0, check4 = 0;

    if (lat > latWest + ((latSouth - latWest)/(lonSouth - lonWest))*(lon-lonWest) )
      check1 = 1;
    if (lat > latSouth + ((latEast - latSouth)/(lonEast - lonSouth))*(lon-lonSouth))
      check2 = 1;
    if (lat < latWest + ((latNorth - latWest)/(lonNorth-lonWest))*(lon-lonWest))
      check3 = 1;
    if (lat < latNorth + ((latEast - latNorth)/(lonEast - lonNorth))*(lon-lonNorth))
      check4 = 1;

    if ((((check1)&&check2)&&check3)&&check4){
      if (num<10000)
        file[num] = (to_string(lat) + ", " + to_string(lon));
      else{
        float prob = 10000.0/(num+1);
        if(distribution1(generator)>prob){
          index = distribution2(generator);
          file[index] = (to_string(lat) + ", " + to_string(lon));
        }
      }
      num++;
    }
  }


  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

  shuffle (file.begin(), file.end(), std::default_random_engine(seed));

  for(int i=0; i<10000; i++)
    f2<<file[i]<<endl;

  f1.close();
  f2.close();

}
#endif
