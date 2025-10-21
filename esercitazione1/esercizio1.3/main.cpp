#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <armadillo>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   int M=100000; //numero di valori da prendere
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

   
///////////////////////////////
//esercizio 1.3   
//////////////////////////////

int Nhit=0;
double d=5.0, l=4.0;

ofstream out ("dati1.3.out");

   for(int i=0; i<N; i++){ //ciclo sui blocchi 

      for(int j=0; j<L; j++){ //ciclo sui passi di ogni blocco 

         double x=rnd.Rannyu(0., d); // estraggo una generica posizione per il centro dell'asta normalizzata su d
         
         //double theta=rnd.Rannyu(0., M_PI);
         
         //estraggo le coordinate dei punti nel quadrato
         double u=rnd.Rannyu(-1.,1.);
         double v=rnd.Rannyu(-1.,1.);

         //verifico che i punti estratti siano nel cerchio
         if(u*u+v*v<1){ //se ricade nel cerchio verifico che l'asta intersechi la linea
               double sin=sqrt(v*v)/sqrt(v*v+u*u);
               if (x+l/(2)*sin>d || x-l/2*sin<0 ){ 
                  Nhit+=1;
               } else {Nhit=Nhit;}
         }else {Nhit=Nhit;
                  j-=1; //se il punto è fuori dal cerchio il lancio non è valido perciò va ripetuto
                  }

      }

         aveL+=2*l*L/(Nhit*d);

        aveL2=pow(aveL,2); //quadrato della media di un blocco (A^2)
        
        aveN+=aveL;
        
        aveN2+=aveL2; //somma dei quadrati delle medie di un blocco
        if (i==0) sigmaN=0;
        else sigmaN=sqrt((aveN2/(i+1)-pow(aveN/(i+1),2))/i);

        out<<i<<"    "<<aveN/(i+1)<<"    "<<sigmaN<<endl;

      
        aveL=0;
        Nhit=0;
        aveL2=0;
        
   }
   out.close();

    cout<<"la media finale e': "<<aveN/100<<endl;
    cout<<"la varianza final e' :"<<sigmaN<<endl;



   return 0;
}

