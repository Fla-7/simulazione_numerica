#ifndef __Chromosome__
#define __Chromosome__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "random.h"
#include "city.h"

using namespace std;
using namespace arma;

class Chromosome {

private:
int n_cities=34; //numero di città 
vector<City> chrom;
//arma:: field<City> chrom(dim);
double p_mutation=0.1; //probabilità di mutare

public: 
Chromosome(bool crcl, Random& rnd); 
Chromosome(int dim) : chrom(dim) { }
Chromosome() { }
bool check();
double get_p_mutation();
int get_size();
City get_city(int i);
void set_city(int i, City city);
double get_length2(); //ritorna L^2
Chromosome swap(Random& rnd); //scmabio di due città vicine per creare la popolazione
//mutazioni
void pair_permutation(Random& rnd); //scambio due città vicine
void shift(Random& rnd); //sposto m città vicine di n posizioni 
void permutation(Random& rnd); //scambio due gruppi di m città 
void inversion(Random& rnd); //inverto l'ordine di m città 
void mutation(Random& rnd, int rank); //assegno le mutazioni a seconda della probabilità 

};



#endif 

