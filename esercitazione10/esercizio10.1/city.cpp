
#include <cmath>
#include <cstdlib>
#include <string>
#include "city.h"

using namespace std;
using namespace arma;

City :: City(bool crcl, Random& rnd ){
  //cout<<"usando costruttore normale di city"<<endl;
  //ofstream out;
  //out.open("coordinate.dat", ios::app);
  if(crcl==true){
    double theta=rnd.Rannyu(0., 2*M_PI); 
    coords.push_back(cos(theta));
    coords.push_back(sin(theta));
    //cout<<coords.at(0)<<"    "<<coords.at(1)<<endl;
  } else{
    coords.push_back(rnd.Rannyu(-1.,1.));
    coords.push_back(rnd.Rannyu(-1.,1.));
  }
//out.close();
}

City::City(double x, double y){
  coords.push_back(x);
  coords.push_back(y);
}

vector <double> City :: get_coords(){
  return coords;
}

/*double City :: get_dist2(vector<double> city1, vector<double> city2){
  return pow(city1[0]-city2[0],2)+pow(city1[1]-city2[1], 2);
}*/

double City :: get_dist2(City& other){
  return pow(coords[0]-other.coords[0],2)+ pow(coords[1]-other.coords[1],2);
}

