#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   int M=100000; //numero di valori da prendere
   int N=100; //numero di blocchi
   int L=M/N; //numero di valori in ogni blocco

   float aveL_call=0, aveL_put=0; //media di un blocco
   float aveL2_call=0, aveL2_put=0; 
   float sumN_call=0, sumN_put=0; //somma delle medie di ogni blocco
   float aveN_call=0, aveN_put=0; //media progressiva
   float aveN2_call=0, aveN2_put=0;
   float err_call=0, err_put=0;

   float S0=100; //prezzo a t=0
   float T=1;  // delivery time
   float K=100;  //strike price
   float r=0.1; //free risk interest rate
   float sigma=0.25; // volatility
   float S=0; //prezzo
   float C_call=0, C_put=0;
   float max_call=0, max_put=0; //per trovare il massimo 
   float Z=0;

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


   rnd.SaveSeed();

   //CALCOLO DIRETTO  

   ofstream outcall ("call_direct.out");
   ofstream output ("put_direct.out");

    for(int i=0; i<N; i++){
      for(int j=0; j<L; j++){
         Z=rnd.Gauss(0,1);
         S=S0*exp((r-pow(sigma,2)/2)*T+sigma*Z*sqrt(T));

         C_call+=exp(-r*T)*fmax(0., S-K);

         C_put+=exp(-r*T)*fmax(0., K-S);
         
      }

      aveL_call=C_call/L; //media del blocco
      aveL_put=C_put/L;

      aveL2_call=pow(aveL_call,2);
      aveL2_put=pow(aveL_put,2);

      sumN_call+=aveL_call;
      sumN_put+=aveL_put;

      aveN_call=sumN_call/(i+1); //media dei blocchi man mano
      aveN_put=sumN_put/(i+1);

      aveN2_call+=aveL2_call;
      aveN2_put+=aveL2_put;

      if(i==0){
         err_call=0;
         err_put=0;
      }else{
         err_call=sqrt((aveN2_call/(i+1)-pow(aveN_call,2))/i);
         err_put=sqrt((aveN2_put/(i+1)-pow(aveN_put,2))/i);
      }

      outcall<<aveN_call<<"\t"<<err_call<<endl;
      output<<aveN_put<<"\t"<<err_put<<endl;

      aveL_call=0; //riporto a zero la media del blocco per rientrare nel ciclo
      aveL_put=0;
      C_call=0;
      C_put=0;

    }

    outcall.close();
    output.close();

    cout<<"call option price: "<<aveN_call<<" +- "<<err_call<<endl;
    cout<<"put option price: "<<aveN_put<<" +- "<<err_put<<endl;

    //azzero le variabili per riusarle nel calcolo successivo
    aveN_call = 0;
    aveN_put = 0;
    aveN2_call = 0;
    aveN2_put = 0;
    sumN_call = 0;
    sumN_put = 0;
    err_call = 0;
    err_put = 0;


   //CALCOLO DISCRETIZZATO 

   outcall.open("call_discrete.out");
   output.open("put_discrete.out");

   for(int i=0; i<N; i++){
      for(int j=0; j<L; j++){
         S=S0;
         for(int k=0; k<100; k++){
            Z=rnd.Gauss(0,1);
            S=S*exp((r-pow(sigma,2)/2)*T/100+sigma*Z*sqrt(T/100));
            
         }

         C_call+=exp(-r*T)*fmax(0., S-K);

         C_put+=exp(-r*T)*fmax(0., K-S);

         
      }

      aveL_call=C_call/L; //media del blocco
      aveL_put=C_put/L;

      aveL2_call=pow(aveL_call,2);
      aveL2_put=pow(aveL_put,2);

      sumN_call+=aveL_call;
      sumN_put+=aveL_put;

      aveN_call=sumN_call/(i+1); //media dei blocchi man mano
      aveN_put=sumN_put/(i+1);

      aveN2_call+=aveL2_call;
      aveN2_put+=aveL2_put;

      //cout<<i<<"\t"<<aveN2_call/(i+1)<<"\t"<<pow(aveN_put,2)<<endl;

      if(i==0){
         err_call=0;
         err_put=0;
      }else{
         err_call=sqrt((aveN2_call/(i+1)-pow(aveN_call,2))/i);
         err_put=sqrt((aveN2_put/(i+1)-pow(aveN_put,2))/i);
      }

      outcall<<aveN_call<<"\t"<<err_call<<endl;
      output<<aveN_put<<"\t"<<err_put<<endl;

      aveL_call=0; //riporto a zero la media del blocco per rientrare nel ciclo
      aveL_put=0;
      C_call=0;
      C_put=0;

    }

    outcall.close();
    output.close();

    cout<<"call option price: "<<aveN_call<<" +- "<<err_call<<endl;
    cout<<"put option price: "<<aveN_put<<" +- "<<err_put<<endl;


   return 0;
}


