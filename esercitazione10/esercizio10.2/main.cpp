#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <armadillo>
#include "random.h"
#include "city.h"
#include "chromosome.h"
#include "population.h"
#include "mpi.h"

using namespace std;
using namespace arma;

int main (int argc, char *argv[]){

   int size, rank; //voglio usare 11 ricerche parallele 
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//leggo le coordinate delle cità
int n_cities=110;
vector<double> x(n_cities);
vector<double>y(n_cities);
ifstream prov("cap_prov_ita.dat");
if(!prov){
   cerr<<"impossibile aprire il file"<<endl;
}else{
   for(int i=0; i<n_cities; i++){
      prov >> x[i]>>y[i];
   }
}
prov.close();

//calcolo la matrice delle distanze (al quadrato per usare L^2!!) per evitare di ripetere il calcolo ogni volta
arma::mat dist_2(n_cities, n_cities);
for(int i=0; i<n_cities; i++){
   for(int k=0; k<n_cities; k++){
      dist_2(i,k)=pow(x.at(i)-x.at(k),2) + pow(y.at(i)-y.at(k),2);
   }
}

//inizializzazione del generatore di numeri casuali con dipendenza dal rank 
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
            seed[0]+=rank; //aggiungo dipendenza dal rank per avere popolazioni diverse
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   bool crcl=false;

   Population pop(n_cities, x, y, rnd); //creo la popolazione totale di tutti i continenti
   
   //cout<<"popolazione creata"<<endl;

   if(pop.check()) //cout<<"popolazione creata con successo"<<endl;
   //cout<<"numero di individui    "<<pop.get_size()<<endl;
   /*for(int i=0; i<pop.get_size(); i++){
      Chromosome chrom=pop.get_element(i);
      City city=chrom.get_city(0);
      vector<double> coords=city.get_coords();
      cout<<i+1<<"      "<<"x:   "<<coords[0]<<"y: "<<coords[1]<<endl;
   }*/


   pop.order(); //ordino gli individui dal migliore al peggiore
   /*for(int i=0; i<pop.get_size(); i++){
      Chromosome chrom=pop.get_element(i);
      cout<<i+1<<"  "<<chrom.get_length()<<endl;
   } // ho verificato che la funzione order funziona*/

   //cout<<"popolazione ordinata"<<endl;

   /*ofstream one("inizio.dat");
   one<<"X:    Y:"<<endl;
   for(int i=0; i<pop.get_element(0).get_size(); i++){
      one<<pop.get_element(0).get_city(i).get_coords().at(0)<<"      "<<pop.get_element(0).get_city(i).get_coords().at(1)<<
            "      "<< pop.get_element(1).get_city(i).get_coords().at(0)<<"      "<<pop.get_element(1).get_city(i).get_coords().at(1)<<endl;
      
   }
one.close();*/


   Population final=pop.evolution(rnd, rank, size);

   //seleziono il percorso migliore tra tutti i rank, gestendo nel rank 0

   struct {
        double length;
        int rank_id;
    } local_best, global_best;
   local_best.length = final.get_element(0).get_length2();
   local_best.rank_id = rank;
   MPI_Allreduce(&local_best, &global_best, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);


   if(rank==0){
      cout << "Migliore lunghezza globale: " << global_best.length << " trovata su Rank " << global_best.rank_id << endl;
   }

   if(rank==global_best.rank_id){
      Chromosome best_chrom;
      best_chrom = final.get_element(0);
      ofstream out("best_path.dat");
      out<<"X:    Y:    "<<global_best.rank_id<<endl;;
      for(int i=0; i<best_chrom.get_size(); i++){
         out<<best_chrom.get_city(i).get_coords().at(0)<<"      "<<best_chrom.get_city(i).get_coords().at(1)<<endl;
      
      }
      out.close();
   
   }
   MPI_Finalize();
   return 0;
}

