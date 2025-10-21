

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "city.h"
#include "chromosome.h"
#include "population.h"

using namespace std;

//funzioni necessarie

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

   bool crcl=true;

   Population pop(crcl, rnd); //creo la popolazione
   cout<<"popolazione creata"<<endl;

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

   cout<<"popolazione ordinata"<<endl;

   ofstream one("inizio.dat");
   one<<"X:    Y:"<<endl;
   for(int i=0; i<pop.get_element(0).get_size(); i++){
      one<<pop.get_element(0).get_city(i).get_coords().at(0)<<"      "<<pop.get_element(0).get_city(i).get_coords().at(1)<<
            "      "<< pop.get_element(1).get_city(i).get_coords().at(0)<<"      "<<pop.get_element(1).get_city(i).get_coords().at(1)<<endl;
      
   }
one.close();
   Population final=pop.evolution(rnd);

   ofstream out("best_path.dat");
   out<<"X:    Y:"<<endl;
   for(int i=0; i<pop.get_element(0).get_size(); i++){
      out<<final.get_element(0).get_city(i).get_coords().at(0)<<"      "<<final.get_element(0).get_city(i).get_coords().at(1)<<endl;
      
   }
   out.close();
////////////////////////////////
//prova  mutazioni
////////////////////////////////

/*pop.get_element(0).pair_permutation(rnd);

ofstream pp("prova_pp.dat");

pp<<"X:    Y:"<<endl;
   for(int i=0; i<pop.get_element(0).get_size(); i++){
      pp<<pop.get_element(0).get_city(i).get_coords().at(0)<<"      "<<pop.get_element(0).get_city(i).get_coords().at(1)<<endl;
      
   }
   pp.close();


pop.get_element(0).shift(rnd);

ofstream shift("prova_shift.dat");

shift<<"X:    Y:"<<endl;
   for(int i=0; i<pop.get_element(0).get_size(); i++){
      shift<<pop.get_element(0).get_city(i).get_coords().at(0)<<"      "<<pop.get_element(0).get_city(i).get_coords().at(1)<<endl;
      
   }
   shift.close();

pop.get_element(0).inversion(rnd);

ofstream inv("prova_inversion.dat");

inv<<"X:    Y:"<<endl;
   for(int i=0; i<pop.get_element(0).get_size(); i++){
      inv<<pop.get_element(0).get_city(i).get_coords().at(0)<<"      "<<pop.get_element(0).get_city(i).get_coords().at(1)<<endl;
      
   }
   inv.close();

   pop.get_element(0).permutation(rnd);

ofstream perm("prova_permutation.dat");

perm<<"X:    Y:"<<endl;
   for(int i=0; i<pop.get_element(0).get_size(); i++){
      perm<<pop.get_element(0).get_city(i).get_coords().at(0)<<"      "<<pop.get_element(0).get_city(i).get_coords().at(1)<<endl;
      
   }
   perm.close();

   Chromosome mom=pop.get_element(0);
   Chromosome dad=pop.get_element(1);

   arma:: field<Chromosome> figli=pop.crossover(mom, dad, rnd); //suggerito da chat ma non ha senso secondo me 
   //figli.at(0)=pop.crossover(mom, dad, rnd).at(0);
   //figli.at(1)=pop.crossover(mom, dad, rnd).at(1);
   City prova=figli.at(0).get_city(0);
cout<<"ciao"<<endl;
   ofstream cross("prova_crossover.dat");

   cross<<"X:    Y:"<<endl;
   for(int i=0; i<pop.get_element(0).get_size(); i++){
      cross<<figli.at(0).get_city(i).get_coords().at(0)<<"      "<<figli.at(0).get_city(i).get_coords().at(1)<<
            "      "<<figli.at(1).get_city(i).get_coords().at(0)<<"      "<<figli.at(1).get_city(i).get_coords().at(1)<<endl;
      
   }
   cross.close();*/


   return 0;
}

