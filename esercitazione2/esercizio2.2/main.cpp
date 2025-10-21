#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   int M=10000; //numero di volte che devo effettuare i RW
   int N=100; //numero di blocchi
   int L=M/N; //numero di valori in ogni blocco
   int P=100; //numero di passi in ogni RW

   double R[3]={0}; //array con le coordinate della posizione 

   double aveL[100]={0}; //media di un blocco
   double aveL2[100]={0}; //quadrato della media di un blocco
   double aveN[100]={0}; //media progressiva sui blocchi
   double aveN2[100]={0}; //media progressiva dei quadrati 
   double sigmaN[100]={0}; //deviazione standard della media 



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


// RW DISCRETO

   ofstream out ("RW_discreto.out");

   for(int i=0; i<N; i++){ //ciclo sui blocchi

      for(int j=0; j<L; j++){ //ciclo dentro al blocco

         for(int k=0; k<P; k++){ //ciclo sui passi del RW
             //scelgo la direzione in cui muovermi tra x,y,z
            int d=round(rnd.Rannyu(0,2)); //arrotondo all'intero più vicino per capire quale coordinata variare 
            if(rnd.Rannyu(-1,1)<0){
               R[d]-=1;
            }else{
               R[d]+=1;
            }
            aveL[k]+=pow(R[0],2)+pow(R[1],2)+pow(R[2],2); //alla fine del ciclo su j, ogni elemento di aveL è la somma di L termini
         }

         //azzero le posizioni per iniziare il RW successivo 
         R[0]=0;
         R[1]=0;
         R[2]=0;

         //qui ho finito un blocco
         
      }

      for(int f=0; f<P; f++){
         aveL[f]/=L;
         aveL[f]=sqrt(aveL[f]); //faccio qui la radice in modo da ottenere l'errore su rad(r2)
         aveL2[f]=pow(aveL[f],2);
         //media e media quadra di ogni blocco vengono accumulati
         aveN[f]+=aveL[f];
         aveN2[f]+=aveL2[f];

         //azzero le medie relative ai singoli blocchi prima di entrare nel blocco successivo
         aveL[f]=0;
         aveL2[f]=0;

      }

       
   //qui ho finito tutti i blocchi 
        
   } 

   //stampa dei risultati

   for(int p=0; p<P; p++){
      sigmaN[p]=sqrt((aveN2[p]/N-pow(aveN[p]/N,2))/(N-1));
      out<<p+1<<"      "<<aveN[p]/N<<"    "<<sigmaN[p]<<endl;

      //azzero gli accumulatori per la simulazione successiva nel continuo
      aveN[p]=0;
      aveN2[p]=0;

   }
 out.close();
 
   
// RW CONTINUO 

   double theta=0;
   double phi=0;

   out.open ("RW_continuo.out");

   for(int i=0; i<N; i++){ //ciclo sui blocchi

      for(int j=0; j<L; j++){ //ciclo dentro al blocco

         for(int k=0; k<P; k++){ //ciclo sui passi del RW
            //estraggo la direzione in cui muovermi
            theta=rnd.Rannyu(0, M_PI);
            phi=rnd.Rannyu(0, 2*M_PI);
            //calcolo le coordinate della posizione dopo lo spostamento
            R[0]+=sin(theta)*cos(phi);
            R[1]+=sin(theta)*sin(phi);
            R[2]+=cos(theta);

            //calcolo la distanza dal centro
            aveL[k]+=pow(R[0],2)+pow(R[1],2)+pow(R[2],2);

         }

         //azzero le posizioni per iniziare il RW successivo 
         R[0]=0;
         R[1]=0;
         R[2]=0;

         //qui ho finito un blocco
         
      }

      for(int f=0; f<P; f++){
         aveL[f]/=L;
         aveL[f]=sqrt(aveL[f]); //faccio qui la radice in modo da ottenere l'errore su rad(r2)
         aveL2[f]=pow(aveL[f],2);
         //media e media quadra di ogni blocco vengono accumulati
         aveN[f]+=aveL[f];
         aveN2[f]+=aveL2[f];

         //azzero le medie relative ai singoli blocchi prima di entrare nel blocco successivo
         aveL[f]=0;
         aveL2[f]=0;

      }

       
   //qui ho finito tutti i blocchi 
        
   } 

   //stampa dei risultati

   for(int p=0; p<P; p++){
      sigmaN[p]=sqrt((aveN2[p]/N-pow(aveN[p]/N,2))/(N-1));
      out<<p+1<<"      "<<aveN[p]/N<<"    "<<sigmaN[p]<<endl;

   }
   
   out.close();

    
   return 0;

}

