#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){


   int M=10000;
   int I[4]={1, 2, 10, 100};
   double S1=0;
   double S2=0;
   double S3=0;

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes.txt");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


   for (int i=0; i<4; i++){
      string filename="distribuzioni_N"+to_string(I[i])+".dat";
     ofstream out (filename);
      for (int j=0; j<M; j++){

         S1=0;
         S2=0;
         S3=0;

         for (int n=1; n<=I[i]; n++){
            S1+=rnd.Rannyu();
            S2+=rnd.Expo(1.);
            S3+=rnd.Lorentz(1., 0.);
         }
         S1=S1/I[i];
         S2=S2/I[i];
         S3=S3/I[i];

         out<<S1<<"     "<<S2<<"    "<<S3<<endl;


      }


      out.close();
   }

   
   


   return 0;
}

