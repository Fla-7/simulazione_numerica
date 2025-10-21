

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

//funzioni necessarie
double err (double sum_ave, double sum_ave2, int i); //calcola l'errore nella media a blocchi
double psi(double x, double mu, double sigma); //valuta la psi in un punto ->mi serve per l'accettazione del metropolis
double Hpsi_psi (double x, double mu, double sigma); //valuta l'azione di H su psi in un punto (quindi la derivata seconda di psi e il potenziale)
bool metro (double x, double y, double mu, double sigma, Random& rnd); //verifica se accettare o no il passo proposto tramite metropolis
double set_delta(double x0, double mu, double sigma, Random& rnd, double n_steps); //ottiene il passo corretto per avere un'accettanza del 50%
//double data_blocking (int n_steps, int n_blk, ); //funzione che stampi le medie progressive in un file e ritorni il valore della media finale

int main (int argc, char *argv[]){

//inizializzazione del generatore di numeri casuali
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

   //variabili per esercizio 8.1
   //double mu=1.; //fisso dei valori per poter calcolare l'integrale, verranno ottimizzati nel 8.2
   //double sigma=0.5;
   double mu=-0.728281;
   double sigma=0.616058;
   double x=0; //punto di partenza
   double y=0.;
   double acceptance=0;
   double n_accepted=0;
   //double delta=0.5;
   double psi_2=0; //per il modulo quadro di psi (esercizio 8.2)

   int n_steps=100000; //numero di passi
   int n_blk=100; //numero di blocchi
   int steps_per_blk=n_steps/n_blk; //numero di passi in ogni blocco
   double ave_blk=0; //media di un blocco
   double ave_prog=0; //media progressiva
   double ave_blk2=0;
   double ave_prog2=0;
   double error=0; //dev std della media 

   double delta=set_delta(x,mu,sigma,rnd, n_steps);
   
   ofstream out ("dati.out");
   out<<"Parametri:  mu="<<mu<<" sigma="<<sigma<<endl;
   out<<"Passo ottenuto per accettazione al 50%:   delta="<<delta<<endl<<endl;
   out<<"#     BLOCK:      H_AVE:       ERROR:     ACCEPTANCE:"<<endl;

   ofstream fz("prob_psi.out");

   for(int i=0; i<n_blk; i++){ //ciclo sui blocchi

      for(int j=0; j<steps_per_blk; j++){ //ciclo dentro al blocco
         y=x+rnd.Rannyu(-1.,1.)*delta;
         if(metro(x,y,mu,sigma,rnd)){
            x=y;
            n_accepted+=1;
         } 
         ave_blk+=Hpsi_psi(x, mu, sigma);
         psi_2=pow(psi(x, mu, sigma),2);
         fz<<x<<"    "<<psi_2<<endl;
      }
      acceptance=(double)n_accepted/(double)steps_per_blk;

      ave_blk=ave_blk/(double)steps_per_blk;
      ave_blk2=pow(ave_blk,2); //quadrato della media di un blocco (A^2)
        
      ave_prog+=ave_blk;
      ave_prog2+=ave_blk2; //somma dei quadrati delle medie di un blocco
      error=err(ave_prog, ave_prog2,i);
      //aveN=aveN/(i+1); //mi serve STAMPARE le medie progressive ma non devo effettivamente modificare la variabile aveN perchè dopo il blocco successivo mi serve che sia solo la somma delle medie dei blocchi, non la media progressiva
      out << setw(12) << i+1
          << setw(12) << ave_prog/double(i+1)
          << setw(12) << error
          << setw(12) << acceptance <<endl;
        
      ave_blk=0; //azzero la media del singolo blocco prima di entrare nel successivo
      ave_blk2=0;
      n_accepted=0;
        
   }

   out.close();

   cout<<"calcolo completato!"<<endl;  

   return 0;
}

double err(double sum_ave, double sum_ave2, int blk){
   double errore;
   if(blk==0){
      errore=0;
   }else{
      errore=sqrt((sum_ave2/(double)(blk+1)-pow(sum_ave/(double)(blk+1),2))/(double)blk);
   }
   return errore;
}

double psi(double x, double mu, double sigma){
   double psi=exp(-pow(x-mu,2)/(2*pow(sigma,2)))+exp(-pow(x+mu,2)/(2*pow(sigma,2)));
   return psi;
}

double Hpsi_psi(double x, double mu, double sigma){
    double sigma2 = sigma * sigma;
    double sigma4 = sigma2 * sigma2;

    double psix = psi(x, mu, sigma);
    if (psix < 1e-12) return 0.0;  // evita problemi numerici

    double d2psi = 
        exp(-pow(x - mu, 2) / (2 * sigma2)) * (pow(x - mu, 2) / sigma4 - 1.0 / sigma2) +
        exp(-pow(x + mu, 2) / (2 * sigma2)) * (pow(x + mu, 2) / sigma4 - 1.0 / sigma2);

    double kinetic = -0.5 * d2psi / psix;
    double potential = pow(x, 4) - 5./2. * pow(x, 2);

    return kinetic + potential;
}


bool metro(double x, double y, double mu, double sigma, Random& rnd){
   bool decision=false;
   double A=0;
   if(1<=pow(psi(y,mu,sigma),2 )/pow(psi(x,mu,sigma),2)) A=1;
   else A=pow(psi(y,mu,sigma),2 )/pow(psi(x,mu,sigma),2);
   double r=rnd.Rannyu();
   if(r<=A) decision=true;
   
   return decision;   
}

double set_delta(double x0, double mu, double sigma, Random& rnd, double n_steps){

    double delta = 1.0;
    int max_attempts = 1000;
    double acceptance = 0;
    double target = 0.5;
    double tolerance = 0.01;
   double x = x0;
    ofstream out("verifying_acceptance.dat");

    for(int attempt = 0; attempt < max_attempts; attempt++){
        int n_accepted = 0;
        
        for(int j = 0; j < n_steps; ++j){
            double y = x + rnd.Rannyu(-1., 1.) * delta;
            if (metro(x, y, mu, sigma, rnd)) {
                x = y;
                n_accepted++;
            }
        }

        acceptance = (double)n_accepted / n_steps;
        out << attempt << " " << delta << " " << acceptance << endl;

        if(fabs(acceptance - target) < tolerance){
            break;
        } else {
            if (acceptance > target) {
                delta *= 1.01;
            } else {
                delta *= 0.99;
            }
        }
    }

    out.close();
    return delta;

   }



