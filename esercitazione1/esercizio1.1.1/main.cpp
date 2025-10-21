
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   int M=10000; //numero di valori da prendere
    int N=100; //numero di blocchi
    int L=M/N; //numero di valori in ogni blocco
    float aveN=0; //media su tutti gli N
    float aveL=0; //media in un blocco
    float aveN2=0;//media dei quadrati
    float aveL2=0; 
   float sigmaN=0; //varianza totale


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


   ofstream out ("dati1.1.1.out");

   for(int i=0; i<N; i++){

      for(int j=0; j<L; j++){
         aveL+=rnd.Rannyu();
         
      }

        aveL=aveL/L;
        aveL2=pow(aveL,2); //quadrato della media di un blocco (A^2)
        
        aveN=aveN+aveL;
        //aveN=aveN/(i+1); //mi serve STAMPARE le medie progressive ma non devo effettivamente modificare la variabile aveN perchè dopo il blocco successivo mi serve che sia solo la somma delle medie dei blocchi, non la media progressiva
        aveN2=(aveN2+aveL2); //somma dei quadrati delle medie di un blocco
        if(i==0) sigmaN=0;
        else sigmaN=sqrt((aveN2/(i+1)-pow(aveN/(i+1),2))/i);

        out<<i<<"    "<<aveN/(i+1)<<"    "<<sigmaN<<endl;

      
        aveL=0;
        aveL2=0;
        
   }

   out.close();
   
    cout<<"la media finale e': "<<aveN/100<<endl;
    cout<<"la varianza final e' :"<<sigmaN<<endl;




   return 0;
}