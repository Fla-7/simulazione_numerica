#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <armadillo>
#include <algorithm>
#include <vector>
#include "random.h"

using namespace std;

//funzioni necessarie
double err (double sum_ave, double sum_ave2, int i); //calcola l'errore nella media a blocchi
double psi(double x, double mu, double sigma); //valuta la psi in un punto ->mi serve per l'accettazione del metropolis
double Hpsi_psi (double x, double mu, double sigma); //valuta l'azione di H su psi in un punto (quindi la derivata seconda di psi e il potenziale)
bool metro (double x, double y, double mu, double sigma, Random& rnd); //verifica se accettare o no il passo proposto tramite metropolis
bool metro_SA(double H_old, double H_new, double T, Random& rnd); //metropolis per l'accettazione dei nuovi parametri
double set_delta(double x0, double mu, double sigma, Random& rnd, double n_steps); //ottiene il passo corretto per avere un'accettanza del 50%
vector<double> exp_H(int n_blk, int n_steps, double T, double mu, double sigma, Random& rnd, bool stampa); //calcola il valore di aspettazione di H e lo stampa in un file

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

   //variabili per esercizio 8.2
   double mu_old=1.; 
   double sigma_old=1.;
   double delta_mu=0.5;
   double delta_sigma=0.5;
   double T_start=2.;
   double T_end=0.01;
   double delta_T=0.99; //diminuisco lentamente la temperatura
   
   int n_steps=10000; //numero di passi
   int n_blk=100; //numero di blocchi
   bool stampa=false;

   vector<double> dati_old=exp_H(n_blk,n_steps, T_start, mu_old,sigma_old,rnd, stampa);//calcolo H con dei parametri mu e sigma
   double H_old=dati_old[0];

   ofstream H_SA("H_SA_steps.dat"); //file con i dati per stampare i valori di H progressivi per ogni step di SA
   ofstream H_opt("optimized_H.dat"); //file con i dati blocco per blocco dell'ultimo passo di SA

   H_SA<<"BETA:        ENERGY:       ERROR:      ACCEPTANCY:"<<endl;
   H_opt.close(); //copio i valori dal file dati.dat -> NON USATO

   for(double T=T_start; T>=T_end; T*=delta_T){ 
      for(int i=0; i<100; i++){ //faccio 100 tentativi per ogni temperatura
      double next_T = T * delta_T;
      bool last_temperature = (next_T < T_end); //se mi trovo nell'ultima temperatura, stampo i valori di H della media a blocchi
      
      //propongo nuovi mu e sigma

      double mu_new=mu_old+rnd.Rannyu(-1.,1.)*T*delta_mu;
      double sigma_new;
            do {
                sigma_new = fabs(sigma_old + rnd.Rannyu(-1, 1) * delta_sigma);
            } while (sigma_new <= 0.1); //evito divisioni per numeri troppo piccoli 
      vector<double> dati_new=exp_H(n_blk, n_steps, T, mu_new, sigma_new, rnd, last_temperature);
      double H_new=dati_new[0];

      //verifico se accettare la mossa
      if(metro_SA(H_old, H_new, T, rnd)){
        mu_old=mu_new;
        sigma_old=sigma_new;
        
   
        H_old=H_new;
      } else{
        dati_old=exp_H(n_blk, n_steps, T, mu_old, sigma_old, rnd, last_temperature);
        H_old=dati_old[0];

      }
      
      
      }
      H_SA<<" "<<T<<"     "<<H_old<<"     "<<dati_old[1]<<"       "<<dati_old[2]<<endl;
   }
   cout<<"calcolo completato!"<<endl;  
   

   return 0;
}


vector<double> exp_H(int n_blk, int n_steps, double T ,double mu, double sigma, Random& rnd, bool stampa){

   //double H=0.;
   vector<double> dati;
   double x=0; //punto di partenza
   double y=0.;
   double delta_H=1.; //passo del metropolis per il calcolo di H
   double acceptance=0;
   double n_accepted=0;

   int steps_per_blk=n_steps/n_blk; //numero di passi in ogni blocco
   double ave_blk=0; //media di un blocco
   double ave_prog=0; //media progressiva
   double ave_blk2=0;
   double ave_prog2=0;
   double error=0; //dev std della media 

   delta_H=set_delta(x,mu,sigma,rnd, n_steps);
   ofstream out ("dati.out");
   if(stampa==true){
      out<<"Parametri: T="<<T<<" mu="<<mu<<" sigma="<<sigma<<endl; 
      out<<"Passo ottenuto per accettazione al 50%:   delta="<<delta_H<<endl<<endl;
      out<<"#     BLOCK:      H_AVE:       ERROR:     ACCEPTANCE:"<<endl;
   }

   int block_eq = 0; // numero di blocchi da scartare per l'equilibrazione


   for(int i=0; i<n_blk; i++){ //ciclo sui blocchi

      for(int j=0; j<steps_per_blk; j++){ //ciclo dentro al blocco
      
         y=x+rnd.Rannyu(-1.,1.)*delta_H;
         if(metro(x,y,mu,sigma,rnd)){
            x=y;
            n_accepted+=1;
         } 
         

         if (i < block_eq)
                continue;
       ave_blk+=Hpsi_psi(x, mu, sigma);
      }
   
      acceptance=(double)n_accepted/(double)steps_per_blk;

      ave_blk=ave_blk/(double)steps_per_blk;
      ave_blk2=pow(ave_blk,2); //quadrato della media di un blocco (A^2)
        
      ave_prog+=ave_blk;
      ave_prog2+=ave_blk2; //somma dei quadrati delle medie di un blocco
      error=err(ave_prog, ave_prog2,i);

      if(stampa==true){
         out << setw(12) << i+1
             << setw(12) << ave_prog/(double)(i+1)
             << setw(12) << error
             << setw(12) << acceptance<<endl;
      }
        
      ave_blk=0; //azzero la media del singolo blocco prima di entrare nel successivo
      ave_blk2=0;
      n_accepted=0;
        
   }

   dati.push_back(ave_prog/(n_blk-block_eq));
   dati.push_back(err(ave_prog, ave_prog2, n_blk-block_eq));
   dati.push_back(acceptance); //verifico l'accettazione con cui è stata calcolata H nell'ultimo blocco

   //return H=ave_prog/n_blk; //stampa come risultato la media finale dopo tutti i blocchi
   return dati;

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

bool metro_SA(double H_old, double H_new, double T, Random& rnd){
   bool decision=false;
   if(H_new<H_old) decision=true;
   else {
        double p=exp(-1./T*(H_new-H_old));
        if(rnd.Rannyu()<p) decision=true;
   };
   
   return decision;   
}

double set_delta(double x0, double mu, double sigma, Random& rnd, double n_steps){
   
    double delta = 1.0;
    int max_attempts = 1000;
    double acceptance = 0;
    double target = 0.5;
    double tolerance = 0.01;

    ofstream out("verifying_acceptance.dat");

    for(int attempt = 0; attempt < max_attempts; ++attempt){
        int n_accepted = 0;
        double x = x0;

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

   